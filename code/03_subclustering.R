library(Seurat)
library(here)
library(ggplot2)
library(harmony)
library(dplyr)
library(patchwork)
library(scCustomize)

plot_dir <- here("plots", "03_subclustering")
processed_dir <- here("processed-data", "03_subclustering")

seurat <- readRDS(here("processed-data","02_seurat_pipeline", "DLS_new_clusters.rds"))
seurat
# An object of class Seurat 
# 43506 features across 28801 samples within 2 assays 
# Active assay: SCT (20574 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap


# ========== Broad cell typing ============
set.seed(125)
# dot plots of broad cell type markers
markers <- c("Snap25", "Syt1",
             "Slc17a7", "Slc17a6",
             "Gad1", "Gad2", "Slc32a1",
             "Sst", "Pvalb",
             "Lamp5", "Vip",
             "Pdyn", "Nts",
             "Npy","Penk",
             "Gfap", "Aqp4",
             "Mbp", "Mobp",
             "Cx3cr1", "P2ry12",
             "Pdgfra","Cspg4",
             "Foxj1","Tmem119")
             

# make blue, white, red color palette
colors_use <- c("blue", "white", "red")


pdf(here(plot_dir, "broad_cell_type_markers.pdf"), width=8, height=8)
Clustered_DotPlot(seurat, features = markers, colors_use_exp = colors_use, cluster_feature=FALSE)
dev.off()

# new broad annotations

broad_annotations <- list(
  "GABA" = c(16,8,14,24,11,4,12,3,5,0,2,1,6,19,10,22),
  "vGlut1" = c(7,15,17,13,23),
  "vGlut2" = c(18,21),
  "Glia" = c(9,20)
)

# Create a vector for broad cell types
seurat$broad_cell_type <- NA  # Initialize with NA

# Assign new labels with sequential numbering for each category
for (cell_type in names(broad_annotations)) {
  clusters <- broad_annotations[[cell_type]]
  counter <- 1  # Initialize counter for numbering
  
  for (cluster in clusters) {
    seurat$broad_cell_type[seurat$seurat_clusters == cluster] <- paste0(cell_type, ".", counter)
    counter <- counter + 1  # Increment counter for the next cluster
  }
}

# Check the distribution
table(seurat$broad_cell_type)
Idents(seurat) <- "broad_cell_type"


# plot with new labels
pdf(here(plot_dir, "broad_cell_type_markers_new_annotations.pdf"), width=8, height=8)
Clustered_DotPlot(seurat, features = markers, colors_use_exp = colors_use, cluster_feature=FALSE)
dev.off()


#  ==== bar plot of broad cell type percentages  ====
broad_cell_type_percentages <- table(seurat$broad_cell_type) / ncol(seurat) * 100
broad_sums <- data.frame(CellType = rownames(broad_sums), Percentage = broad_sums$broad_sums, row.names = NULL)

# Save plot to PDF
pdf(here(plot_dir, "broad_cell_type_percentages.pdf"), width=2.5, height=3.5)
ggplot(broad_sums, aes(x = CellType, y = Percentage, fill=CellType)) +
  geom_bar(stat = "identity") +
  labs(title = "Broad Cell Type Percentages",
       x = "Broad Cell Type",
       y = "Percentage") +
  ggpubr::theme_pubr() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.x=element_blank())
dev.off()


# ==== stacked bars of percent samples that make up each broad cell type ====
sample_percents <- table(seurat$broad_cell_type, seurat$Sample_ID)
#          DLS1 DLS2 DLS3
# GABA.1    198  182  188
# GABA.10  1446 1421 1427
# GABA.11   837  710  802
# GABA.12  1128  993 1061
# GABA.13   535  502  480

# Convert counts to row-wise percentages
sample_percents <- table(seurat$broad_cell_type, seurat$Sample_ID)
sample_percents_norm <- prop.table(sample_percents, margin = 1) * 100  # Convert to percentage

# Convert matrix to dataframe while keeping row and column names
df_long <- as.data.frame(sample_percents_norm) %>%
  tibble::rownames_to_column() 
# rowname     Var1 Var2     Freq
# 1        1   GABA.1 DLS1 34.85915
# 2        2  GABA.10 DLS1 33.67490
# 3        3  GABA.11 DLS1 35.63218
# 4        4  GABA.12 DLS1 35.44940
# 5        5  GABA.13 DLS1 35.26697
# 6        6  GABA.14 DLS1 33.60825

# rename Var1 to CellType and Var2 to Sample
colnames(df_long) <- c("", "Cell_Type", "Sample", "Percentage")

# Create stacked bar plot
pdf(here(plot_dir, "broad_cell_type_sample_percentages.pdf"), width=9, height=3)
ggplot(df_long, aes(x = Cell_Type, y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Percentage of Samples in Each Broad Cluster",
       x = "Broad Cluster",
       y = "Percentage",
       fill = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
dev.off()


# save final broad cell type clusters
saveRDS(seurat, here(processed_dir, "DLS_broad_celltypes.rds"))


# ========= Drop vGlut1 neurons and re-clustering ======

# Ensure broad_cell_type is a character vector
seurat$broad_cell_type <- as.character(seurat$broad_cell_type)

# Identify cells that are NOT vGlut1
cells_to_keep <- !grepl("^vGlut1", seurat$broad_cell_type)

# Subset the Seurat object to keep only these cells
seurat <- seurat[, cells_to_keep]

# Verify that vGlut1 neurons are removed
table(seurat$broad_cell_type)
# GABA.1  GABA.10  GABA.11  GABA.12  GABA.13  GABA.14  GABA.15  GABA.16   GABA.2   GABA.3   GABA.4   GABA.5   GABA.6 
# 568     4294     2349     3182     1517      485     1134      301     1202      592       67      929     1693 
# GABA.7   GABA.8   GABA.9   Glia.1   Glia.2 vGlut2.1 vGlut2.2 
# 891     2057     1618     1140      426      504      420 

# === Rerun seurat pipeline with SCTransform ===
seurat <- SCTransform(seurat, vars.to.regress = "percent.mito", verbose = FALSE)

seurat <- RunPCA(seurat, verbose = FALSE) %>%
  RunHarmony("Sample_ID", plot_convergence=TRUE)

seurat <- FindNeighbors(seurat, reduction="harmony", dims=1:30)
seurat <- FindClusters(seurat, resolution=0.5)
seurat
# An object of class Seurat 
# 43293 features across 25369 samples within 2 assays 
# Active assay: SCT (20361 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap

Idents(seurat) <- "seurat_clusters"

#  ======== UMAP =========
seurat <- RunUMAP(seurat, reduction="harmony", dims=1:20)
seurat

# plot umap
pdf(here(plot_dir, "new_umap_clusters.pdf"), width=5, height=5)
DimPlot(seurat, reduction="umap", label=TRUE)
dev.off()

# umap with samples
pdf(here(plot_dir, "new_umap_samples.pdf"), width=5, height=5)
DimPlot(seurat, reduction="umap", group.by="Sample_ID")
dev.off()

# umap with markers
pdf(here(plot_dir, "new_umap_markers.pdf"), width=8, height=8)
FeaturePlot(seurat, reduction="umap", features=c("Snap25", "Mbp", "Slc17a7", 'Gad1'))
dev.off()

# umap with qc metrics
pdf(here(plot_dir, "new_umap_qc.pdf"), width=12, height=4)
FeaturePlot(seurat, reduction="umap", features=c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol=3)
dev.off()



# ========== Broad cell typing ============

# dot plots of broad cell type markers
markers <- c("Snap25", "Syt1",
             "Slc17a7", "Slc17a6",
             "Gad1", "Gad2", "Slc32a1",
             "Sst", "Pvalb",
             "Lamp5", "Vip",
             "Pdyn", "Nts",
             "Npy","Penk",
             "Gfap", "Aqp4",
             "Mbp", "Mobp",
             "Cx3cr1", "P2ry12",
             "Pdgfra","Cspg4",
             "Foxj1","Tmem119")


# make blue, white, red color palette
colors_use <- c("blue", "white", "red")


pdf(here(plot_dir, "new_broad_cell_type_markers.pdf"), width=8, height=8)
Clustered_DotPlot(seurat, features = markers, colors_use_exp = colors_use, cluster_feature=FALSE)
dev.off()


# new broad annotations
broad_annotations <- list(
  "GABA" = c(16,14,22,23,1,0,20,12,5,2,7,10,9,13,3,6,4,19),
  "vGlut2" = c(17,15),
  "Glia" = c(11,8,21,18)
)

# Create a vector for broad cell types
seurat$broad_cell_type <- NA  # Initialize with NA

# Assign new labels with sequential numbering for each category
for (cell_type in names(broad_annotations)) {
  clusters <- broad_annotations[[cell_type]]
  counter <- 1  # Initialize counter for numbering
  
  for (cluster in clusters) {
    seurat$broad_cell_type[seurat$seurat_clusters == cluster] <- paste0(cell_type, ".", counter)
    counter <- counter + 1  # Increment counter for the next cluster
  }
}

# Check the distribution
table(seurat$broad_cell_type)
Idents(seurat) <- "broad_cell_type"


# plot with new labels
pdf(here(plot_dir, "new_broad_cell_type_markers_with_annotations.pdf"), width=8, height=8)
Clustered_DotPlot(seurat, features = markers, colors_use_exp = colors_use, cluster_feature=FALSE)
dev.off()

table(seurat$broad_cell_type)


# save
saveRDS(seurat, here(processed_dir, "DLS_new_clusters_without_vGlut1.rds"))

# ====== GABA annotations =======

# subset to only GABA neurons
seurat.gaba <- seurat[, grepl("^GABA", seurat$broad_cell_type)]
seurat.gaba
# An object of class Seurat 
# 43293 features across 21657 samples within 2 assays 
# Active assay: SCT (20361 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap


# ====== subclustering GABA neurons ======
set.seed(123)
seurat.gaba <- SCTransform(seurat.gaba, vars.to.regress = "percent.mito", verbose = FALSE)

seurat.gaba <- RunPCA(seurat.gaba, verbose = FALSE) %>%
  RunHarmony("Sample_ID", plot_convergence=TRUE)

seurat.gaba <- FindNeighbors(seurat.gaba, reduction="harmony", dims=1:20)
seurat.gaba <- FindClusters(seurat.gaba, resolution=0.5, algorithm=4)
seurat.gaba
# An object of class Seurat 
# 43293 features across 25369 samples within 2 assays 
# Active assay: SCT (20361 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap

Idents(seurat.gaba) <- "seurat_clusters"

#  ======== UMAP =========
seurat.gaba <- RunUMAP(seurat.gaba, reduction="harmony", dims=1:20)
seurat.gaba

# plot umap
pdf(here(plot_dir, "Gaba_umap_subclusters.pdf"), width=5, height=5)
DimPlot(seurat.gaba, reduction="umap", label=TRUE)
dev.off()


# === Ploting ===
gaba_markers <- c("Sst", "Pvalb",  "Lamp5","Pdyn", "Nts","Penk", "Prkcd")

# make blue, white, red color palette
colors_use <- c("blue", "white", "red")


pdf(here(plot_dir, "GABA_subclusters_cell_type_markers.pdf"), width=8, height=5)
Clustered_DotPlot(seurat.gaba, features = gaba_markers, colors_use_exp = colors_use, cluster_feature=FALSE)
dev.off()



# umap with qc metrics
pdf(here(plot_dir, "Gaba_subcluster_clusters.pdf"), width=5, height=5)
DimPlot(seurat.gaba, reduction="umap", label=TRUE)
dev.off()

# umap with qc metrics
pdf(here(plot_dir, "Gaba_subcluster_umap_Sst.pdf"), width=5, height=5)
FeaturePlot(seurat.gaba, reduction="umap", features="Sst")
dev.off()

pdf(here(plot_dir, "Gaba_subcluster_umap_Nts.pdf"), width=5, height=5)
FeaturePlot(seurat.gaba, reduction="umap", features="Nts")
dev.off()

pdf(here(plot_dir, "Gaba_subcluster_umap_Npy.pdf"), width=5, height=5)
FeaturePlot(seurat.gaba, reduction="umap", features="Npy")
dev.off()

pdf(here(plot_dir, "Gaba_subcluster_umap_Penk.pdf"), width=5, height=5)
FeaturePlot(seurat.gaba, reduction="umap", features="Penk")
dev.off()

pdf(here(plot_dir, "Gaba_subcluster_umap_Pdyn.pdf"), width=5, height=5)
FeaturePlot(seurat.gaba, reduction="umap", features="Pdyn")
dev.off()

pdf(here(plot_dir, "Gaba_subcluster_umap_Trpc4.pdf"), width=5, height=5)
FeaturePlot(seurat.gaba, reduction="umap", features="Trpc4")
dev.off()


# save
saveRDS(seurat.gaba, here(processed_dir, "DLS_gaba_subclusters.rds"))
                        
