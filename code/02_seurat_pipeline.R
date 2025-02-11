library(Seurat)
library(here)
library(ggplot2)
library(harmony)
library(dplyr)
library(patchwork)

plot_dir <- here("plots", "02_seurat_pipeline")
processed_dir <- here("processed-data", "02_seurat_pipeline")

load(here("processed-data","01_quality_control", "DLSclusters_qc.Robj"))
seurat
# An object of class Seurat 
# 22932 features across 28801 samples within 1 assay 
# Active assay: RNA (22932 features, 2356 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: ica, tsne

# ========= SCTransform ==========
seurat <- SCTransform(seurat, vars.to.regress = "percent.mito", verbose = FALSE)
# An object of class Seurat 
# 43506 features across 28801 samples within 2 assays 
# Active assay: SCT (20574 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: ica, tsne


# ========= Dim Reduction + Harmony Integration ==========
seurat <- RunPCA(seurat, verbose = FALSE) %>%
  RunHarmony("Sample_ID", plot_convergence=TRUE)
seurat
# An object of class Seurat 
# 43506 features across 28801 samples within 2 assays 
# Active assay: SCT (20574 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 5 dimensional reductions calculated: ica, tsne, pca, umap, harmony

# plot PCA before and after intergration
pdf(here(plot_dir, "pca_before_after_harmony.pdf"), width=8, height=4)
p1 <- DimPlot(seurat, reduction="pca", group.by="Sample_ID") + ggtitle("PCA before Harmony")
p2 <- DimPlot(seurat, reduction="harmony", group.by="Sample_ID") + ggtitle("PCA after Harmony")
p1+p2
dev.off()

# harmony pca colored by qc metrics
pdf(here(plot_dir, "pca_harmony_qc.pdf"), width=12, height=4)
FeaturePlot(seurat, reduction="harmony", features=c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol=3)
dev.off()

pdf(here(plot_dir, "pca_harmony_markers.pdf"), width=8, height=8)
FeaturePlot(seurat, reduction="harmony", features=c("Snap25", "Mbp", "Slc17a7", 'Gad1'))
dev.off()

pdf(here(plot_dir, "pca_eblow.pdf"), width=5, height=5)
ElbowPlot(seurat)
dev.off()

# ======== Clustering ==========
seurat@reductions$ica <- NULL
seurat@reductions$tsne <- NULL
seurat@reductions$umap <- NULL

seurat <- FindNeighbors(seurat, reduction="harmony", dims=1:30)
seurat <- FindClusters(seurat, resolution=0.5)
seurat
# An object of class Seurat 
# 43506 features across 28801 samples within 2 assays 
# Active assay: SCT (20574 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, harmony

Idents(seurat) <- "seurat_clusters"

#  ======== UMAP =========
seurat <- RunUMAP(seurat, reduction="harmony", dims=1:20)
seurat
# An object of class Seurat 
# 43506 features across 28801 samples within 2 assays 
# Active assay: SCT (20574 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap

# plot umap
pdf(here(plot_dir, "umap_clusters.pdf"), width=5, height=5)
DimPlot(seurat, reduction="umap", label=TRUE)
dev.off()

# umap with samples
pdf(here(plot_dir, "umap_samples.pdf"), width=5, height=5)
DimPlot(seurat, reduction="umap", group.by="Sample_ID")
dev.off()

# umap with markers
pdf(here(plot_dir, "umap_markers.pdf"), width=8, height=8)
FeaturePlot(seurat, reduction="umap", features=c("Snap25", "Mbp", "Slc17a7", 'Gad1'))
dev.off()

# umap with qc metrics
pdf(here(plot_dir, "umap_qc.pdf"), width=12, height=4)
FeaturePlot(seurat, reduction="umap", features=c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol=3)
dev.off()

# save
saveRDS(seurat, here(processed_dir, "DLS_new_clusters.Robj"))


