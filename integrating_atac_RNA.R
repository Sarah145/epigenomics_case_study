library(Seurat)
library(Signac)
library(tidyverse)
library(patchwork)

atac <- readRDS('tum_norm_integrated.Rds')
rna <- readRDS('tumour_norm_scRNA_labelled.Rds')
p1 <- DimPlot(rna, group.by = "cell_type", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
p1 + p2


# quantify gene activity
gene.activities <- GeneActivity(atac, features = VariableFeatures(rna))

# add gene activities as a new assay
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)
atac <- ScaleData(atac, features = rownames(atac))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$cell_type,
                                     weight.reduction = atac[["integrated_lsi"]], dims = 2:30)

atac <- AddMetaData(atac, metadata = celltype.predictions)

DimPlot(atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
DimPlot(atac, group.by = "sampletype", label = F) 

atac <- FindNeighbors(atac, reduction = "integrated_lsi", dims = 2:30)
atac <- FindClusters(object = atac, graph = 'peaks_snn', verbose = FALSE, algorithm = 3, resolution = 0.5)

ggplot(atac@meta.data, aes(x = predicted.id, fill = peaks_snn_res.0.5)) +
  geom_bar(position = 'fill') +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggplot(atac@meta.data, aes(x = peaks_snn_res.0.5, fill = predicted.id)) +
  geom_bar(position = 'fill') +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

saveRDS(atac, 'tumour_norm_atac_annotated.Rds')
