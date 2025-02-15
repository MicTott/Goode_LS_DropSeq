library(Seurat)
library(here)
library(ggplot2)
library(harmony)
library(dplyr)
library(patchwork)
library(scCustomize)

plot_dir <- here("plots", "04_marker_genes")
processed_dir <- here("processed-data", "04_marker_genes")

seurat <- readRDS(here("processed-data","03_subclustering", "DLS_gaba_subclusters.rds"))
seurat

# ========== Marker genes ============
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
ls.gaba.markers <- FindAllMarkers(seurat, only.pos = TRUE)
top10_markers <- ls.gaba.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()


pdf(here(plot_dir, "Heatmap_marker_genes.pdf"), width=10, height=10)
DoHeatmap(seurat, features = top10_markers$gene, group.by = "seurat_clusters") + NoLegend()
dev.off()

top50_markers <- ls.gaba.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 50) %>%
  ungroup()

# export to csv
write.csv(top50_markers, here(processed_dir, "top50_markers_per_GABAclass.csv"))

# ======= Cell type proportions ========

# Extract metadata
meta_df <- seurat@meta.data

# Calculate cluster proportions per sample
cluster_props <- meta_df %>%
  group_by(Sample_ID, seurat_clusters) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Sample_ID) %>%
  mutate(Proportion = Count / sum(Count))

# Get mean proportions for ordering and labeling
cluster_means <- cluster_props %>%
  group_by(seurat_clusters) %>%
  summarise(Mean_Prop = mean(Proportion) * 100, .groups = "drop") %>%
  arrange(desc(Mean_Prop))

# Order clusters by descending mean proportion
cluster_props$seurat_clusters <- factor(cluster_props$seurat_clusters, levels = cluster_means$seurat_clusters)

# Plot bar plot with mean % text at the top
pdf(here(plot_dir, "Cell_type_proportions.pdf"), width=8, height=3)
ggplot(cluster_means, aes(x = seurat_clusters, y = Mean_Prop, fill=seurat_clusters)) +
  geom_bar(stat = "identity") +  # Bar plot
  geom_text(aes(label = sprintf("%.1f%%", Mean_Prop)), vjust = -0.5, size = 3) +  # Mean % text
  theme_minimal() +
  labs(x = "GABAergic Cell Type Cluster", y = "Mean Proportion (%)", title = "GABAergic Cell Type Proportions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggpubr::theme_pubr() +
  theme(legend.position = "none") +
  ylim(0, 15)
dev.off()


# ======= Denrogram =======
seurat <- BuildClusterTree(seurat)
PlotClusterTree(seurat)


# Visualize co-expression of two features simultaneously
FeaturePlot(seurat, features = c("Sst", "Pdyn"), blend = TRUE)



# Define genes of interest
genes_of_interest <- c("Sst", "Pdyn", "Penk")

# Extract expression data for selected genes
expr_matrix <- GetAssayData(seurat, slot = "data")[genes_of_interest, ]

# Convert to data frame and add cluster labels
expr_df <- as.data.frame(t(as.matrix(expr_matrix)))
expr_df$Cluster <- seurat$seurat_clusters

# Convert expression to binary (1 if expressed, 0 otherwise)
expr_df[, genes_of_interest] <- (expr_df[, genes_of_interest] > 0) * 1  # Convert to binary

# Function to compute % co-expression matrix for each cluster
get_coexpression_percentage <- function(subset_df, genes) {
  total_cells <- nrow(subset_df)  # Total cells in cluster
  
  # Create a data frame to store pairwise co-expression percentages
  coexpression_matrix <- expand.grid(Gene1 = genes, Gene2 = genes, stringsAsFactors = FALSE)
  coexpression_matrix$Coexpression <- NA
  
  # Compute % of cells co-expressing each pair of genes
  for (i in seq_len(nrow(coexpression_matrix))) {
    g1 <- coexpression_matrix$Gene1[i]
    g2 <- coexpression_matrix$Gene2[i]
    
    coexpressed_cells <- sum(subset_df[[g1]] == 1 & subset_df[[g2]] == 1)  # Count cells expressing both
    coexpression_matrix$Coexpression[i] <- (coexpressed_cells / total_cells) * 100  # Convert to percentage
  }
  
  return(coexpression_matrix)
}

# Apply function to each cluster
coexpression_data <- expr_df %>%
  group_by(Cluster) %>%
  group_split() %>%
  lapply(function(df) {
    cluster_coexp <- get_coexpression_percentage(df, genes_of_interest)
    cluster_coexp$Cluster <- unique(df$Cluster)  # Assign cluster ID
    return(cluster_coexp)
  }) %>%
  bind_rows()  # Combine results into a single dataframe

# Convert cluster and genes to factors for ordering
coexpression_data$Cluster <- factor(coexpression_data$Cluster, levels = unique(coexpression_data$Cluster))
coexpression_data$Gene1 <- factor(coexpression_data$Gene1, levels = genes_of_interest)
coexpression_data$Gene2 <- factor(coexpression_data$Gene2, levels = genes_of_interest)

# Add diagonal values with light grey
diagonal_data <- data.frame(
  Gene1 = genes_of_interest,
  Gene2 = genes_of_interest,
  Coexpression = -1,  # Dummy value for diagonal
  Cluster = rep(unique(coexpression_data$Cluster), each = length(genes_of_interest))
)

# Combine with main dataset
coexpression_data <- bind_rows(coexpression_data, diagonal_data)

# Cap co-expression values at 5% for visualization
coexpression_data$Coexpression <- pmin(coexpression_data$Coexpression, 5)

# Plot the heatmap
pdf(here(plot_dir, "Coexpression_heatmap.pdf"), width=10, height=10)
ggplot(coexpression_data, aes(x = Gene1, y = Gene2, fill = Coexpression)) +
  geom_tile(color = "black") +  # Add black borders to tiles
  scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 5), na.value = "lightgrey") +  # Cap at 5% and set diagonal to grey
  facet_wrap(~Cluster, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(size = 12, face = "bold"))
dev.off()