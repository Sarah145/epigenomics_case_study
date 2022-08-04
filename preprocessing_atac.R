library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(tidyverse)
set.seed(1234)

# metadata file contains barcodes for cells passing QC
# downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE181062&format=file&file=GSE181062%5FRCC%5FeightPt%5Fcombined%5FscATAC%2ECD45%5Fplus%5Fcell%5FMetaData%2Etxt%2Egz
tum_norm_md <- read.csv('GSE181062_RCC_eightPt_combined_scATAC.CD45_plus_cell_MetaData.txt.gz', sep = '\t') %>% 
  filter(str_detect(barcode, '(pt1002300_RCC)|(pt1002300_adjN)'))
tum_norm_cells <- paste0(ifelse(tum_norm_md$sampletype == 'RCC', 'RCC-', 'Norm-'), str_extract(tum_norm_md$barcode, '[^#]+$'))


# files extracted from GSE181062_RAW, downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE181062&format=file
# tumour and normal sample fragments were concatenated using following commands
# gunzip GSM5483811_1279543_fragments.tsv.gz GSM5483809_1251334_fragments.tsv.gz
# awk -i inplace '$4="Norm-"$4' GSM5483811_1279543_fragments.tsv
# awk -i inplace '$4="RCC-"$4' GSM5483809_1251334_fragments.tsv
# cat GSM5483809_1251334_fragments.tsv GSM5483811_1279543_fragments.tsv > tumour_normal_fragments.tsv
# sed -i 's/\s/\t/g' tumour_normal_fragments.tsv
# bedtools sort -i tumour_normal_fragments.tsv > tumour_normal_fragments_sort.tsv
# bgzip tumour_normal_fragments_sort.tsv
# tabix -s 1 -b 2 -e 3 tumour_normal_fragments_sort.tsv.gz


# read in fragments
tum_norm_fragments <- CreateFragmentObject(
  path = 'GSE181062_RAW/tumour_normal_fragments_sort.tsv.gz',
  validate.fragments = FALSE,
  cells = tum_norm_cells
)

# call peaks with MACS3
tum_norm_peaks <- CallPeaks(tum_norm_fragments, 
                            macs2.path = '/home/sarah/anaconda3/bin/macs3')

# filter out peaks in blacklist regions
# blacklist came from here: http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/
blacklist <- read.csv('hg38.blacklist.bed.gz', sep='\t', header = F)
for(i in 1:nrow(blacklist)){
  tum_norm_peaks <- tum_norm_peaks[!(seqnames(tum_norm_peaks) == blacklist$V1[i] & tum_norm_peaks@ranges@start >= blacklist$V2[i] & (tum_norm_peaks@ranges@start + tum_norm_peaks@ranges@width) <= blacklist$V3[i]), ]
}

# count fragments in peaks
tum_norm_counts <- FeatureMatrix(
  tum_norm_fragments,
  tum_norm_peaks,
  cells = tum_norm_cells,
  sep = c("-", "-"),
  verbose = TRUE
)

tum_norm_chrom_assay <- CreateChromatinAssay(
  counts = tum_norm_counts,
  sep = c("-", "-"),
  genome = 'hg38',
  fragments = 'GSE181062_RAW/tumour_normal_fragments_sort.tsv.gz',
  min.cells = 10,
  min.features = 200
)

rownames(tum_norm_md) <- tum_norm_cells
tum_norm <- CreateSeuratObject(
  counts = tum_norm_chrom_assay,
  assay = "peaks",
  meta.data = tum_norm_md)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(tum_norm) <- annotations

tum_norm <- RunTFIDF(tum_norm)
tum_norm <- FindTopFeatures(tum_norm, min.cutoff = 'q0')
tum_norm <- RunSVD(tum_norm)

DepthCor(tum_norm) DepthCor(tum_norm) # component 1 correlates with depth

tum_norm <- RunUMAP(object = tum_norm, reduction = 'lsi', dims = 2:30)
tum_norm <- FindNeighbors(object = tum_norm, reduction = 'lsi', dims = 2:30)
tum_norm <- FindClusters(object = tum_norm, verbose = FALSE, algorithm = 3, resolution = 0.5)
DimPlot(object = tum_norm, label = TRUE) + NoLegend()
DimPlot(object = tum_norm, label = TRUE, group.by = 'sampletype') # can see batch effect between samples here
gene.activities <- GeneActivity(tum_norm)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
tum_norm[['RNA']] <- CreateAssayObject(counts = gene.activities)
tum_norm <- NormalizeData(
  object = tum_norm,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(tum_norm$nCount_RNA)
)

# integrate samples
integration.anchors <- FindIntegrationAnchors(
  object.list = SplitObject(tum_norm, split.by = "sampletype"),
  anchor.features = rownames(tum_norm),
  reduction = "rlsi",
  dims = 2:30,
  assay = rep('peaks', 2)
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = tum_norm[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
DimPlot(integrated, group.by = "sampletype") # batch effect looks better

# cluster and look at activity of marker genes to label clusters
integrated <- FindNeighbors(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated <- FindClusters(integrated, algorithm = 4, resolution = 0.2, graph.name = 'peaks_snn')
DimPlot(integrated, group.by = c('peaks_snn_res.0.2'), label = T, repel = T) + NoLegend()
DefaultAssay(integrated) <- 'RNA'
FeaturePlot(
  object = integrated,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

cluster_names <- c('1 - NK cells (tumour + normal)', '2 - T cells (tumour)', '3 - T cells (tumour + normal)', '4 - Phagocytic cells', '5 - Cytotoxic cells', '6 - Unknown', '7 - NK cells (tumour)', '8 - Unknown (normal)', '9 - B cells', '10 - Tregs (tumour)')
integrated$cluster <- cluster_names[as.numeric(integrated$peaks_snn_res.0.2)]
integrated$cluster <- factor(integrated$cluster, levels = cluster_names)

saveRDS(integrated, 'tum_norm_integrated_atac.Rds')
