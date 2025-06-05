library(Seurat)
library(here)
library(ggplot2)
library(harmony)
library(dplyr)
library(patchwork)
library(scCustomize)
library(Nebulosa)

plot_dir <- here("plots", "06_sst_subclustering")
processed_dir <- here("processed-data", "06_sst_subclustering")

seurat <- readRDS(here("processed-data","03_subclustering", "DLS_gaba_subclusters.rds"))
seurat
# An object of class Seurat 
# 43506 features across 28801 samples within 2 assays 
# Active assay: SCT (20574 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap



# subset to SST clusters
clusters <- c(3,4,7,8,16,17,20)
obj.sst <- subset(seurat, idents = clusters)
obj.sst
# An object of class Seurat 
# 42976 features across 6496 samples within 2 assays 
# Active assay: SCT (20044 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap

# ====== Subclustering ======
set.seed(631)
obj.sst <- SCTransform(obj.sst, vars.to.regress = "percent.mito", verbose = FALSE)

obj.sst <- RunPCA(obj.sst, verbose = FALSE) %>%
  RunHarmony("Sample_ID", plot_convergence=TRUE)

obj.sst <- FindNeighbors(obj.sst, reduction="harmony", dims=1:30)
obj.sst <- FindClusters(obj.sst, resolution=0.3, dims = 1:40, k.param=40)
obj.sst
# An object of class Seurat 
# 41041 features across 6496 samples within 2 assays 
# Active assay: SCT (18109 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap


Idents(obj.sst) <- "seurat_clusters"

#  ======== UMAP =========
obj.sst <- RunUMAP(obj.sst, reduction="harmony", dims=1:40)
obj.sst

# plot umap
pdf(here(plot_dir, "UMAP_SST_subclusters.pdf"), width=5, height=5)
DimPlot(obj.sst, reduction="umap", label=TRUE)
dev.off()


# collapse clusters 0 and 5
obj.sst$seurat_clusters <- as.character(obj.sst$seurat_clusters)
obj.sst$seurat_clusters[obj.sst$seurat_clusters == "5"] <- "0"
obj.sst$seurat_clusters[obj.sst$seurat_clusters == "3"] <- "4"
obj.sst$seurat_clusters[obj.sst$seurat_clusters == "1"] <- "9"

# numbers now skip 5. Reassign numbers 0-9
clusters <- as.character(obj.sst$seurat_clusters)

# get the old IDs (e.g. "0","1","2","3","4","6","7","8","9","10")
old_ids <- sort(unique(clusters))

# create the new IDs 0:(n-1) — here 0:9
new_ids <- as.character(0:(length(old_ids)-1))

# remap into a factor with levels 0–9
obj.sst$seurat_clusters <- factor(clusters,
                                  levels = old_ids,
                                  labels = new_ids)

levels(obj.sst$seurat_clusters)
# [1] "0" "1" "2" "3" "4" "5" "6" "7" "8" "9"


Idents(obj.sst) <- "seurat_clusters"

# umap
pdf(here(plot_dir, "UMAP_SST_subclusters_collapsed.pdf"), width=5, height=5)
DimPlot(obj.sst, reduction="umap", label=TRUE) +
  scale_color_manual(values = scCustomize::ColorBlind_Pal()) +
  theme_void()
dev.off()


table(obj.sst$seurat_clusters)
# 0    1    2    3    4    5    6    7 
# 1818  151  836 1447  372  362  336 1174 


# ====== Marker genes =======
ls.gaba.markers <- FindAllMarkers(obj.sst, only.pos = TRUE)
top50_markers <- ls.gaba.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 50) %>%
  ungroup()

# export top 50 markers per cell type
write.csv(top50_markers, here(processed_dir, "SST_subcluster_markers.csv"), row.names = FALSE)



top10_markers <- ls.gaba.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()


pdf(here(plot_dir, "Heatmap_marker_genes.pdf"), width=10, height=10)
DoHeatmap(obj.sst, features = top10_markers$gene, group.by = "seurat_clusters", group.colors=scCustomize::ColorBlind_Pal()) +
  NoLegend() +
  scale_fill_gradientn(colors = c("white", "lightgrey", "red"))
dev.off()


pdf(here(plot_dir, "UMAPs_SSTsubcluster_wNebulosa.pdf"), width=20, height=4)
plots <- plot_density(obj.sst, c("Sst","Nts", "Npy", "Penk", "Pdyn"), raster=FALSE)

# Ensure theme modification is applied to each individual plot
plots <- lapply(plots, function(p) {
  p + theme_void() +
    # increase title size
    theme(plot.title = element_text(size = 20)) 
})

wrap_plots(plots, nrow = 1) # Ensure it uses patchwork's wrap_plots to arrange

dev.off()



# save thing but a normal express umap (DimPlot)
pdf(here(plot_dir, "UMAPs_SSTsubcluster_dimplot.pdf"), width=20, height=4)
FeaturePlot(obj.sst, features = c("Sst","Nts", "Npy", "Penk", "Pdyn"), cols = c("lightgrey", "red"), order = TRUE, ncol=5)
dev.off()







# dotplot of SST subcluster marker genes 
features <- c("Cntn5", "Fam19a1", "Kcnt2", "Car10", "Cpa6", "Gpc5", "Sox6","Nfib")
pdf(here(plot_dir, "Dotplot_SST_subcluster_markers.pdf"), width=6, height=3.25)
DotPlot(obj.sst, features ,cols = c("lightgrey", "red"),) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() +
 # revesrse y order
  scale_y_discrete(limits = rev) 
dev.off()




# ======== UMAP and Marker genes for all gaba dataset =========


# save thing but a normal express umap (DimPlot)
pdf(here(plot_dir, "UMAPs_all_Neurons_dimplot.pdf"), width=20, height=4)
FeaturePlot(seurat, features = c("Sst","Nts", "Npy", "Penk", "Pdyn"), cols = c("lightgrey", "red"), order = TRUE, ncol=5)
dev.off()



# ====== Marker genes =======
ls.gaba.markers <- FindAllMarkers(seurat, only.pos = TRUE)
top10_markers <- ls.gaba.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()


pdf(here(plot_dir, "Heatmap_marker_genes_all.pdf"), width=10, height=10)
DoHeatmap(seurat, features = top10_markers$gene, group.by = "seurat_clusters", group.colors=scCustomize::ColorBlind_Pal()) +
  NoLegend() +
  scale_fill_gradientn(colors = c("white", "lightgrey", "red"))
dev.off()


# ===== Ried et al heatmap ======
library(ComplexHeatmap)
library(circlize)

# make orig assay active
DefaultAssay(seurat) <- "RNA"

genes <- c("Trpc4", "Meis2", "Foxp2",  "Zeb2","Homer2", 
           "Hunk",
           "Sst", "Nts", "Tacr3",
           "Rorb",
           "Tmem132c", "Pax6",
           "Moxd1", "L3mbtl4", "Cux2",
           "Ndnf", 
           "Slc5a7", "Hopx",
           "Cartpt", "Igf2bp2",
           "Esr1", "Calcr",
           "Emilin2",
           "Chst8",
           "Glp1r", 
           "Etv1",
           "Tshz2","Pth2r",
           "Dach2", "Lhx2","Crhr2",
           "Drd3",
           "Met", 
           "Lamp5",
           "Pcsk6"
            )
           
           
# drop "g" form idents ("g1", "g2", etc)
seurat$seurat_clusters <- as.character(seurat$seurat_clusters)
seurat$seurat_clusters <- gsub("g", "", seurat$seurat_clusters)

Idents(seurat) <- "seurat_clusters"

# heatmap using complex heatmap

# 3. Make sure the data are scaled for those features
#    (this will write into seurat[["RNA"]]@scale.data)
seurat <- ScaleData(seurat, features = genes, verbose = FALSE)

# 4. Compute the average scaled expression per identity class
#    (assumes your identities are set to what you want on the heatmap columns)
avg_list <- AverageExpression(
  seurat,
  assays = "RNA",
  slot   = "scale.data",
  features = genes,
  verbose = FALSE
)
avg_scaled <- avg_list$RNA

avg_scaled <- scale(t(avg_scaled))

# drop "g" from row names
rownames(avg_scaled) <- gsub("g", "", rownames(avg_scaled))

# 6. Define a diverging color map centered at 0
col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

# make row annotations with cluster colors
cluster_labels <- sort(as.numeric(unique(Idents(seurat)))-1)
cluster_colors <-  DiscretePalette_scCustomize(num_colors = 22, palette = "polychrome") # Assign colors
names(cluster_colors) <- unique(cluster_labels)

row_ha <- rowAnnotation(
  Cluster = rownames(avg_scaled),
  col = list(Cluster = cluster_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE
)


# 7. Draw the heatmap
p <- Heatmap(
  avg_scaled,
  name                   = "Scaled\nExpr",
  col                    = col_fun,
  cluster_rows           = TRUE,
  cluster_columns        = TRUE,
  show_row_names         = TRUE,
  show_column_names      = TRUE,
  row_names_gp           = gpar(fontsize = 10),
  column_names_gp        = gpar(fontsize = 10),
  heatmap_legend_param   = list(
    title = "Z-score"
  ),
  right_annotation= row_ha
)

pdf(here(plot_dir, "Heatmap_Ried_et_al.pdf"), width=8, height=5)
draw(p, heatmap_legend_side = "right")
dev.off()



# ===== Dot plot instead =====
# re roder cluster levels 
seurat$seurat_clusters <- factor(seurat$seurat_clusters)

new_levels <- levels(seurat$seurat_clusters) %>%
  as.character() %>% 
  as.numeric() %>% 
  sort() %>% 
  as.character()

seurat$seurat_clusters <- factor(
  seurat$seurat_clusters,
  levels = new_levels
)

Idents(seurat) <- seurat$seurat_clusters

# dot plot instead
pdf(here(plot_dir, "Dotplot_Ried_et_al.pdf"), width=10, height=5)
DotPlot(seurat, features = genes, cols = c("lightgrey", "red"), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_gradientn(colors = c("white", "lightgrey", "red"))
dev.off()