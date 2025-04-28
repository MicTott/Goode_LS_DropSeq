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


pdf(here(plot_dir, "dotplot_broad_markers.pdf"), width=8, height=8)
Clustered_DotPlot(seurat, features = markers, colors_use_exp = colors_use, cluster_feature=FALSE)
dev.off()


# $broad_labels with just those 4 broad annotations to the seurat object
seurat$broad_celltype <- NA
seurat$broad_celltype[seurat$seurat_clusters %in% c(16,8,14,24,11,4,12,3,5,0,2,1,6,19,10,22)] <- "GABA"
seurat$broad_celltype[seurat$seurat_clusters %in% c(7,15,17,13,23)] <- "vGlut1"
seurat$broad_celltype[seurat$seurat_clusters %in% c(18, 21)] <- "vGlut2"
seurat$broad_celltype[seurat$seurat_clusters %in% c(9,20)] <- "Glia"

# Check the broad cell type amounts
table(seurat$broad_celltype)
# GABA   Glia vGlut1 vGlut2 
# 22879   1566   3432    924 

Idents(seurat) <- "broad_celltype"


# umap with broad cell type labels
pdf(here(plot_dir, "umap_broad_celltypes.pdf"), width=5, height=5)
DimPlot(seurat, reduction="umap", group.by="broad_celltype")
dev.off()



#  ==== bar plot of broad cell type percentages  ====
broad_sums <- table(seurat$broad_celltype)
broad_celltype_percentages <- table(seurat$broad_celltype) / ncol(seurat) * 100
broad_sums <- data.frame(CellType = rownames(broad_sums), Percentage = broad_sums, row.names = NULL)

# Save plot to PDF
pdf(here(plot_dir, "broad_celltype_percentages.pdf"), width=2.5, height=3.5)
ggplot(broad_sums, aes(x = CellType, y = Percentage.Freq, fill=CellType)) +
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
sample_percents <- table(seurat$broad_celltype, seurat$Sample_ID)
#          DLS1 DLS2 DLS3
# GABA.1    198  182  188
# GABA.10  1446 1421 1427
# GABA.11   837  710  802
# GABA.12  1128  993 1061
# GABA.13   535  502  480

# Convert counts to row-wise percentages
sample_percents <- table(seurat$broad_celltype, seurat$Sample_ID)
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
pdf(here(plot_dir, "broad_celltype_sample_percentages.pdf"), width=9, height=3)
ggplot(df_long, aes(x = Cell_Type, y = Percentage, fill = Cell_Type)) +
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
set.seed(1382)
# Ensure broad_celltype is a character vector
seurat$broad_celltype <- as.character(seurat$broad_celltype)

# Identify cells that are NOT vGlut1
cells_to_keep <- !grepl("^vGlut1", seurat$broad_celltype)

# Subset the Seurat object to keep only these cells
seurat <- seurat[, cells_to_keep]

# Verify that vGlut1 neurons are removed
table(seurat$broad_celltype)
# GABA   Glia vGlut2 
# 22879   1566    924 

# === Rerun seurat pipeline with SCTransform ===
options(future.globals.maxSize = 1000 * 1024^2)
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
pdf(here(plot_dir, "umap_clusters_noGlut1.pdf"), width=5, height=5)
DimPlot(seurat, reduction="umap", label=TRUE)
dev.off()

# umap with samples
pdf(here(plot_dir, "umap_samples_noGlut1.pdf"), width=5, height=5)
DimPlot(seurat, reduction="umap", group.by="Sample_ID")
dev.off()

# umap with broad cell type
pdf(here(plot_dir, "umap_broad_celltypes_noGlut1.pdf"), width=5, height=5)
DimPlot(seurat, reduction="umap", group.by="broad_celltype")
dev.off()

# umap with markers
pdf(here(plot_dir, "umap_markers_noGlut1.pdf"), width=8, height=8)
FeaturePlot(seurat, reduction="umap", features=c("Snap25", "Mbp", "Slc17a6", 'Gad1'))
dev.off()

# umap with qc metrics
pdf(here(plot_dir, "new_umap_qc.pdf"), width=12, height=4)
FeaturePlot(seurat, reduction="umap", features=c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol=3)
dev.off()


# save
saveRDS(seurat, here(processed_dir, "DLS_broad_clusters_without_vGlut1.rds"))






# ====== GABA annotations =======

# subset to only GABA neurons
seurat.gaba <- seurat[, grepl("^GABA", seurat$broad_celltype)]
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

seurat.gaba <- FindNeighbors(seurat.gaba, reduction="harmony", dims=1:40)
seurat.gaba <- FindClusters(seurat.gaba, resolution=0.5, dims = 1:40, k.param=30)
seurat.gaba
# An object of class Seurat 
# 43293 features across 25369 samples within 2 assays 
# Active assay: SCT (20361 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap

Idents(seurat.gaba) <- "seurat_clusters"

#  ======== UMAP =========
seurat.gaba <- RunUMAP(seurat.gaba, reduction="harmony", dims=1:40)
seurat.gaba

# plot umap
pdf(here(plot_dir, "Gaba_umap_subclusters.pdf"), width=5, height=5)
DimPlot(seurat.gaba, reduction="umap", label=TRUE)
dev.off()


# === Ploting ===
gaba_markers <- c("Sst", "Pvalb",  "Lamp5","Pdyn", "Nts","Penk", "Prkcd", "Vip")

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






# ==== stacked bars of percent samples that make up each broad cell type ====
sample_percents <- table(seurat.gaba$seurat_clusters, seurat.gaba$Sample_ID)


# Convert counts to row-wise percentages
sample_percents <- table(seurat.gaba$seurat_clusters, seurat.gaba$Sample_ID)
sample_percents_norm <- prop.table(sample_percents, margin = 1) * 100  # Convert to percentage

# Convert matrix to dataframe while keeping row and column names
df_long <- as.data.frame(sample_percents_norm) %>%
  tibble::rownames_to_column() 


# rename Var1 to CellType and Var2 to Sample
colnames(df_long) <- c("", "Cell_Type", "Sample", "Percentage")

# Create stacked bar plot
pdf(here(plot_dir, "broad_celltype_sample_percentages.pdf"), width=7, height=3)
ggplot(df_long, aes(x = Cell_Type, y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Percentage of Samples in Each Broad Cluster",
       x = "Broad Cluster",
       y = "Percentage",
       fill = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
dev.off()







# save
saveRDS(seurat.gaba, here(processed_dir, "DLS_gaba_subclusters.rds"))
                        
