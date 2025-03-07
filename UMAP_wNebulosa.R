library(Seurat)
library(here)
library(ggplot2)
library(harmony)
library(dplyr)
library(patchwork)
library(scCustomize)
library(circlize)
library(ComplexHeatmap)
library(SingleCellExperiment)
library(Nebulosa)

plot_dir <- here("plots", "05_nebulosa_UMAPs")
processed_dir <- here("processed-data", "04_nebulosa_UMAPs")

seurat <- readRDS(here("processed-data","03_subclustering", "DLS_gaba_subclusters.rds"))
seurat

# convert to SCE
sce <- as.SingleCellExperiment(seurat)
sce
# class: SingleCellExperiment 
# dim: 20044 21657 
# metadata(0):
#   assays(2): counts logcounts
# rownames(20044): 0610007P14Rik 0610009B22Rik ... Timm8a2 Vmn2r73
# rowData names(0):
#   colnames(21657): DLS1_ACCCACTCAACAACCT-1 DLS1_TAGAGCTCATGGTCTA-1 ... DLS3_CACTCCAGTTAGTGGG-1
# DLS3_CGTCAGGGTCTTGCGG-1
# colData names(10): nCount_RNA nFeature_RNA ... broad_cell_type ident
# reducedDimNames(3): PCA HARMONY UMAP
# mainExpName: SCT
# altExpNames(1): RNA


# ======= Plotting with Nebulosa ======

pdf(here(plot_dir, "UMAPs_wNebulosa.pdf"), width=20, height=4)
plots <- plot_density(seurat, c("Sst","Nts", "Npy", "Penk", "Pdyn"), raster=FALSE)

# Ensure theme modification is applied to each individual plot
plots <- lapply(plots, function(p) {
  p + theme_void() +
    # increase title size
    theme(plot.title = element_text(size = 20)) 
})

wrap_plots(plots, nrow = 1) # Ensure it uses patchwork's wrap_plots to arrange

dev.off()


# Penk and Nts
pdf(here(plot_dir, "UMAPs_wNebulosa_Sst_Penk.pdf"), width=4, height=4)
plots<-plot_density(seurat, c("Sst", "Penk"), raster=FALSE, joint=TRUE)

# Ensure theme modification is applied to each individual plot
plots <- lapply(plots, function(p) {
  p + theme_void() +
    # increase title size
    theme(plot.title = element_text(size = 20)) 
})

# only plot te last plot
wrap_plots(plots[[length(plots)]], nrow = 1) # Ensure it uses patchwork's wrap_plots to arrange
dev.off()


# Penk and Pdyn
pdf(here(plot_dir, "UMAPs_wNebulosa_Penk_Pdyn.pdf"), width=4, height=4)
plots<-plot_density(seurat, c("Penk", "Pdyn"), raster=FALSE, joint=TRUE)

# Ensure theme modification is applied to each individual plot
plots <- lapply(plots, function(p) {
  p + theme_void() +
    # increase title size
    theme(plot.title = element_text(size = 20)) 
})

wrap_plots(plots[[length(plots)]], nrow = 1) # Ensure it uses patchwork's wrap_plots to arrange
dev.off()


# Sst and Npy
pdf(here(plot_dir, "UMAPs_wNebulosa_Sst_Npy.pdf"), width=4, height=4)
plots<-plot_density(seurat, c("Sst", "Npy"), raster=FALSE, joint=TRUE)

# Ensure theme modification is applied to each individual plot
plots <- lapply(plots, function(p) {
  p + theme_void() +
    # increase title size
    theme(plot.title = element_text(size = 20)) 
})

wrap_plots(plots[[length(plots)]], nrow = 1) # Ensure it uses patchwork's wrap_plots to arrange
dev.off()



# subset to these clusters and make a heatmap of these genes
# SST+/PDYN+: 1, 8
# PDYN+/PENK+: 14
# PENK+: 7
# SST+/NTS+: 11
# NTS+:12
# NTS+/PENK: 5
# SST+/NPY+: 15

# subset to these clusters
clusters <- c(1, 2, 8, 14, 7, 11, 12, 5, 15, 17)
sce.subset <- sce[, sce$ident %in% clusters]



library(SingleCellExperiment)
library(ComplexHeatmap)
library(circlize)
library(scuttle)  # For aggregateAcrossCells
library(here)

# Define marker genes
genes <- c("Sst", "Nts", "Npy", "Penk", "Pdyn")

# Aggregate expression per cluster using scuttle
sce.agg <- aggregateAcrossCells(sce.subset, 
                                ids = colData(sce.subset)$ident,
                                statistics="sum",
                                use.assay.type="logcounts")

# Extract log-normalized counts of selected genes from aggregated object
expr_grouped <- logcounts(sce.agg)[genes, ]

# Min-Max Scaling per gene
min_max_scale <- function(x) {
  if (max(x) == min(x)) return(rep(0, length(x)))  # Handle zero-variance cases
  (x - min(x)) / (max(x) - min(x))
}
expr_scaled <- apply(expr_grouped, 1, min_max_scale)
expr_scaled <- as.matrix(expr_scaled)  # Ensure matrix format



# ---- Calculate the percentage of nuclei per cluster ----
cluster_counts <- table(colData(sce)$ident)  # Count cells per cluster
cluster_percentage <- (cluster_counts / sum(cluster_counts)) * 100  # Convert to percentage

# get only clusters we need
cluster_percentage <- as.matrix(cluster_percentage[clusters])
# [,1]
# 1  12.753382
# 8   5.102276
# 14  2.747380
# 7   6.215081
# 11  3.952533
# 12  3.887888
# 5   7.775777
# 15  2.627326


# reorder the rows of cluster percentage so they are in numerical order
cluster_percentage <- cluster_percentage[order(as.numeric(rownames(cluster_percentage))),]


# ---- Create a Barplot Annotation ----
bar_anno <- HeatmapAnnotation(
  GABA_Percentage = anno_barplot(cluster_percentage, 
                                   gp = gpar(fill = "gray40"), 
                                   height = unit(1, "cm"))
)



# Save heatmap as PNG
pdf(here(plot_dir, "heatmap_marker_genes.pdf"), width=5.5, height=3.5)

Heatmap(t(expr_scaled),
        name = "Scaled Expression",
        col = colorRamp2(c(0, 1), c("white", "red")),  # White to red scale
        cluster_rows = TRUE, 
        cluster_columns = FALSE,  # Do not cluster clusters
        show_row_names = TRUE,
        show_column_names = TRUE,
        top_annotation = bar_anno,  # Add barplot annotation,
        column_title = "Cluster",
        column_title_side = "bottom")

dev.off()
