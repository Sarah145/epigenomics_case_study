library(Seurat)
library(tidyverse)
library(Matrix)

# files downloaded from https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/query/acc.cgi?acc=GSE181061
mtx <- readMM('GSE181061_ccRCC_4pt_scRNAseq_CD45plus_matrix.mtx.gz')
barcodes <- read.csv('GSE181061_ccRCC_4pt_scRNAseq_CD45plus_barcodes.tsv.gz', sep='\t', header = F)
md <- read.csv('GSE181061_ccRCC_4pt_scRNAseq_CD45plus_final_Metadata.txt.gz', sep = '\t')
genes <- read.csv('GSE181061_ccRCC_4pt_scRNAseq_CD45plus_genes.tsv.gz', sep='\t', header = F)
rownames(mtx) <- genes$V1
colnames(mtx) <- barcodes$V1

pt_tn_barcodes <- md %>% filter(Patient == '1002300', tissue != 'pbmc') %>% pull(barcode)
pt_tn_mtx <- mtx[,pt_tn_barcodes]
pt_tn_md <- md %>% filter(Patient == '1002300', tissue != 'pbmc')
rownames(pt_tn_md) <- pt_tn_md$barcode

seu <- CreateSeuratObject(counts = pt_tn_mtx,
                          meta.data = pt_tn_md,
                          min.cells = 10,
                          min.features = 200)
seu$rep <- ifelse(str_detect(seu$barcode, '_1_'), 1, 2)
seu$sample <- paste0(seu$tissue, '_', seu$rep)
seu.list <- SplitObject(seu, split.by = "sample")

# normalize and identify variable features for each dataset independently
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seu.list)
seu_anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = features)
seu.combined <- IntegrateData(anchorset = seu_anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seu.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
seu.combined <- ScaleData(seu.combined, verbose = FALSE)
seu.combined <- RunPCA(seu.combined, npcs = 15, verbose = FALSE)
seu.combined <- RunUMAP(seu.combined, reduction = "pca", dims = 1:15)
seu.combined <- FindNeighbors(seu.combined, reduction = "pca", dims = 1:15)
seu.combined <- FindClusters(seu.combined, resolution = 0.5)
# Visualization
DimPlot(seu.combined, reduction = "umap", group.by = "sample")
DimPlot(seu.combined, reduction = "umap", group.by = "tissue")
DimPlot(seu.combined, reduction = "umap", label = T)
FeaturePlot(seu.combined,
            features = c('MS4A1', 'CD8A', 'NKG7', 'LYZ', 'S100A8', 'HBA1'),
            ncol = 2)

degs <- FindAllMarkers(seu.combined,
                       assay = 'RNA',
                       test.use = 'LR',
                       latent.vars = 'sample')

markers <- degs %>% group_by(cluster) %>% filter(p_val_adj <= 0.05) %>%
  filter(avg_log2FC > 0) %>%
  slice_min(order_by = p_val_adj, n = 20) %>%
  slice_max(order_by = avg_log2FC, n = 5)

ggplot(markers, aes(x = gene, y = -log10(p_val_adj +1))) +
  geom_text(aes(label = gene), angle = 90) +
  facet_wrap(~cluster, nrow = 5, scales = 'free')

DimPlot(seu.combined, group.by = 'seurat_clusters', label = T)

# b cells
FeaturePlot(seu.combined,
            features = c('MS4A1', 'CD79A'))
VlnPlot(seu.combined,
        feature = c('MS4A1', 'CD79A'), ncol = 2)
# # clusters 7 and 11 are B cells - will see if I can further delineate them
# b_degs <- FindMarkers(seu.combined, assay = 'RNA', ident.1 = '7', ident.2 = '11')
# FeaturePlot(seu.combined,
#             features = c('TCL1A', 'IL4R'))
# I think cluster 7 may have more T cell contamination than cluster 11 but both are B cells
cell_types <- rep('x', ncol(seu.combined))
cell_types[seu.combined$seurat_clusters == '9'] <- 'B cells'

# monocytes and DCs
FeaturePlot(seu.combined,
            features = c('CST3', 'LYZ', 'S100A8', 'FCER1A'))
FeaturePlot(seu.combined,
            features = c('CD14', 'FCGR3A')) # cd14 vs cd16 monocytes
FeaturePlot(seu.combined,
            features = c('CLEC10A', 'IL3RA', 'IRF8')) # cDC vs pDC
cell_types[seu.combined$seurat_clusters == '7'] <- 'CD14 monocytes'
cell_types[seu.combined$seurat_clusters == '10'] <- 'CD16 monocytes'
cell_types[seu.combined$seurat_clusters == '12'] <- 'Dendritic cells'

# T cells
FeaturePlot(seu.combined,
            features = c('CD3E', 'CD8A', 'CD4', 'CCR7'))
FeaturePlot(seu.combined,
            features = c('IL7R', 'CCL5', 'CTLA4'))

cell_types[seu.combined$seurat_clusters == '5' | seu.combined$seurat_clusters == '3'] <- 'CD8 memory T cells' 
cell_types[seu.combined$seurat_clusters == '1'] <- 'CD8 naive T cells' 
cell_types[seu.combined$seurat_clusters == '4'] <- 'CD4 T cells' 
cell_types[seu.combined$seurat_clusters == '14'] <- 'Regulatory T cells' 

# NK cells
FeaturePlot(seu.combined,
            features = c('NKG7', 'KLRB1', 'GNLY', 'PRF1'))
cell_types[seu.combined$seurat_clusters %in% c('0', '6', '11', '13')] <- 'NK cells'
cell_types[seu.combined$seurat_clusters == '2'] <- 'NKT cells' # expresses T and NK cell markers

# Leftover cell types
cell_types[seu.combined$seurat_clusters == '8'] <- 'Dying cells'
cell_types[seu.combined$seurat_clusters %in% c('15', '16')] <- 'Unknown'

sum(cell_types %in% c('Dying cells', 'Unknown'))
seu.combined$cell_type <- cell_types
Idents(seu.combined) <- seu.combined$cell_type

DimPlot(seu.combined, label = T)

# remove dying and unknown cells
seu.combined <- seu.combined[,seu.combined$cell_type %ni% c('Unknown', 'Dying cells')]
DimPlot(seu.combined, label = T)
saveRDS(seu.combined, 'tumour_norm_scRNA_labelled.Rds')
