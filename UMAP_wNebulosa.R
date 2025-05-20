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
# SST+/PDYN+: 4, 8
# PDYN+/PENK+: 13
# PENK+: 12
# SST+/NTS+: 20
# NTS+:3
# NTS+/PENK: 3
# SST+/NPY+: 17
# SST+: 16, 7

# subset to these clusters
clusters <- c(4,8,13,12,20,3,17,7, 16)
sce.subset <- sce[, sce$ident %in% clusters]



library(SingleCellExperiment)
library(ComplexHeatmap)
library(circlize)
library(scuttle)  # For aggregateAcrossCells
library(here)

# Define marker genes
genes <- c("Sst", "Nts", "Npy", "Penk", "Pdyn", "Glp1r", "Crhr2", "Pnoc", "Mchr1")

# Aggregate expression per cluster using scuttle
sce.agg <- aggregateAcrossCells(sce.subset, 
                                ids = colData(sce.subset)$ident,
                                statistics="sum",
                                use.assay.type="logcounts")


# drop unused idents
sce.agg$ident <- droplevels(sce.agg$ident)

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
cluster_percentage <- as.matrix(cluster_percentage[clusters+1])
# [,1]
# 4  6.7276169
# 8  5.1900078
# 13 2.7381447
# 12 3.4769359
# 20 0.6879993
# 3  8.6438565
# 12 3.4769359
# 17 1.5514614
# 15 2.1194071
# 7  5.4624371


# reorder the rows of cluster percentage so they are in numerical order
cluster_percentage_ordered <- cluster_percentage[order(as.numeric(rownames(cluster_percentage))),]


# ---- Create a Barplot Annotation ----
bar_anno <- HeatmapAnnotation(
  GABA_Percentage = anno_barplot(cluster_percentage_ordered, 
                                   gp = gpar(fill = "gray40"), 
                                   height = unit(1, "cm"))
)



# Save heatmap as PNG
pdf(here(plot_dir, "heatmap_marker_genes_subset_otherGenes.pdf"), width=5.5, height=3.5)

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




# ========== Find top 10 markers for these clusters =============

# Find markers for each cluster using scran
library(scran)

# drop unused levels
sce.subset$ident <- droplevels(sce.subset$ident)
unique(sce.subset$ident)

markers <- findMarkers(sce.subset, 
                              groups = sce.subset$ident, 
                              test.type = "wilcox", 
                              direction = "up")

markers[[1]][1:10]
# DataFrame with 20044 rows and 10 columns
# Top      p.value          FDR summary.AUC     AUC.4     AUC.7     AUC.8    AUC.12    AUC.13    AUC.16
# <integer>    <numeric>    <numeric>   <numeric> <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
#   Auts2          1 2.16885e-210 5.43405e-207    0.809827  0.809827  0.749443  0.623046  0.828651  0.556926  0.642300
# Celf2          1 1.41476e-167 2.18134e-164    0.972531  0.703581  0.586031  0.475693  0.812375  0.430790  0.670024
# Cntn5          1 1.62742e-293 3.26200e-289    0.882083  0.778889  0.882083  0.860396  0.387982  0.433829  0.640043
# Dlgap1         1 3.03798e-200 6.08933e-197    0.910457  0.552520  0.703151  0.661604  0.646224  0.910457  0.636461

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)    # you only need map()

top10_df <- 
  # 1) turn your SimpleList into a tibble of list-columns
  as.list(markers) %>%                  
  enframe(name = "cell_type", value = "df") %>%  
  
  # 2) grab the top ten, flatten the rownames into a column
  mutate(
    df = map(df, ~ head(.x, 10) %>% 
               as.data.frame(stringsAsFactors = FALSE) %>%
               rownames_to_column("marker"))
  ) %>%
  
  # 3) row-bind them all, auto-filling any missing columns with NA
  unnest(cols = df)

# Voilà:
top10_df

library(scater)

# 1) Pull out your top-10 marker genes
marker_genes <- unique(top10_df$marker)

# 2) Plot the heatmap of mean expression per cluster
plotGroupedHeatmap(
  sce.subset,
  features = marker_genes,   # the genes to show
  group  = "ident",
  scale=TRUE,
  center=TRUE,
  cluster_rows=FALSE,
  cluster_cols=FALSE
  # which colData column to average over
)


library(scater)
library(scuttle)
library(ComplexHeatmap)
library(circlize)        # for colorRamp2

# (a) subset to markers
sce.sub <- sce.subset[marker_genes, ]

# (b) pseudobulk (mean here)
sce.pseudo <- aggregateAcrossCells(
  sce.sub,
  ids        = sce.sub$ident,
  statistics = "mean"
)

# (c) extract counts and subset
expr_mat <- assay(sce.pseudo, "counts")[marker_genes, , drop = FALSE]

# (d) min–max scale rows to [0,1]
min_max_scale <- function(x) (x - min(x)) / (max(x) - min(x))
scaled_mat <- t(apply(expr_mat, 1, min_max_scale))
scaled_mat <- scaled_mat[marker_genes, ]

# (e) transpose so clusters → rows, markers → columns
mat_t <- t(scaled_mat)

# (e) same Heatmap call as above
col_fun <- colorRamp2(c(0, 1), c("white", "red"))

# (g) draw
pdf(here(plot_dir, "heatmap_top10_markers.pdf"), width=10, height=4)
Heatmap(
  mat_t,
  col = col_fun,
  name            = "Scaled expr.",
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  row_title       = "Cluster",
  column_title    = "Markers"
)
dev.off()
