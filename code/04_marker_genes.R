library(Seurat)
library(here)
library(ggplot2)
library(harmony)
library(dplyr)
library(patchwork)
library(scCustomize)
library(circlize)
library(ComplexHeatmap)

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



selected_genes <- c()  # Track genes that have been assigned

top1_markers <- ls.gaba.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(cluster) %>%  # Ensure original order
  group_split() %>%  # Process each cluster separately
  purrr::map_df(~ {
    df <- .x %>% dplyr::filter(!gene %in% selected_genes)  # Remove previously assigned genes
    if (nrow(df) > 0) {
      selected_genes <<- c(selected_genes, df$gene[1])  # Store chosen gene
      return(df[1, ])  # Keep only the first available gene per cluster
    } else {
      return(NULL)  # If no unique gene remains for this cluster, it is skipped
    }
  })

markers <- top1_markers$gene


pdf(here(plot_dir, "Violins_TopMarker.pdf"), width=5, height=8)
VlnPlot_scCustom(seurat, features = markers, pt.size = 0.1,  stack=TRUE, flip = TRUE) + NoLegend() 
dev.off()




# ======= Heatmap of top 5 marker genes per cell type =======

top5_markers <- ls.gaba.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup()

markers <- top5_markers$gene

# Step 2: Aggregate expression per cluster
# Calculate average expression per cluster for the marker genes
avg_expr <- AverageExpression(seurat, features = markers, return.seurat = FALSE)
expr_matrix <- avg_expr$RNA  # Extract matrix from list

# drop g from colnames
colnames(expr_matrix) <- gsub("g", "", colnames(expr_matrix))

# Step 3: Scale the data for visualization
scaled_expr <- t(scale(t(expr_matrix)))  # Row-wise scaling

# Step 4: Define cluster annotation
cluster_labels <- factor(colnames(scaled_expr))  # Cluster names
cluster_colors <-  DiscretePalette_scCustomize(num_colors = 22, palette = "polychrome") # Assign colors
names(cluster_colors) <- unique(cluster_labels)



# Step 5: Generate the heatmap
pdf(here(plot_dir, "Heatmap_top5_marker_genes.pdf"), width=6, height=9)
Heatmap(scaled_expr,
        name = "Expression",
        cluster_rows = FALSE,  # Cluster genes
        cluster_columns = FALSE,  # Don't cluster cell type clusters
        show_column_names = TRUE,
        column_title = "Cell Type Clusters",
        row_title = "Top 5 Marker Genes",
        #col = colorRamp2::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),  # Color scale
        top_annotation = HeatmapAnnotation(Cluster = cluster_labels, 
                                           col = list(Cluster = cluster_colors))
)
dev.off()



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
# change clusters names from"g1", "g2", ... to "1", "2", etc, ...
seurat$seurat_clusters <- gsub("g", "", seurat$seurat_clusters)

# Set the cleaned cluster labels as factor
seurat$seurat_clusters <- factor(seurat$seurat_clusters, levels = sort(unique(seurat$seurat_clusters)))
Idents(seurat) <- seurat$seurat_clusters
seurat@tools$BuildClusterTree <- NULL

# Rebuild the tree using cleaned identities
seurat <- BuildClusterTree(seurat, dims = 1:20)  # adjust dims as needed
PlotClusterTree(seurat)



# === clustered dotplots of dev markers ===
library(scCustomize)

# Define genes of interest
markers <- c("Trpc4", "Homer2", "Foxp2", "Meis2", "Zeb2")

# make a white to red color scale
my_colors <- colorRamp2(c(0, 1), c("white", "red"))

# Create a custom dot plot
pdf(here(plot_dir, "DotPlot_developmentMarkers.pdf"), width=9, height=4)
Clustered_DotPlot(seurat, features = markers, colors_use_exp=viridis::rocket(n = 20, direction = -1))
dev.off()


# ======== UMAP with better colors ========

pdf(here(plot_dir, "UMAP_with_scCustomize_colors.pdf"), width=5, height=5)
DimPlot_scCustom(seurat) +
  theme_void()
dev.off()








 # =========== Co expression stuff ===========
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




# ========= Dot plots of marker genes from LS paper ==========
# https://www.biorxiv.org/content/10.1101/2024.02.15.580381v1.full.pdf

genes <- c("Foxp2", "Gucy1a1", "Drd3", "Rorb", "Tshz2", "Chst8", "Pth2r", "Zeb2", 
           "Lamp5", "Glp1r", "Ndnf", "Esr1", "Cux2", "Tacr3", "Calcr", "Met", 
           "Emilin2", "Pcsk6", "Igf2bp2", "Pax6", "Tmem132c", "Dach2", "Lhx2", 
           "Crhr2", "Hunk", "Nts", "Sst", "Cartpt", "Hopx", "Etv1", "Chat", 
           "Moxd1", "Slc5a7", "L3mbtl4", "Meis2", "Sfta3-ps")

pdf(here(plot_dir, "DotPlot_marker_genes_Reid_et_al.pdf"), width=6, height=8)
DotPlot_scCustom(seurat, features = genes,  flip_axes=TRUE, colors_use = c("blue", "white", "red")) + NoLegend()
dev.off()



genes <- c("Meis2","Foxp2", "Zeb2", 
           "Nts", "Sst","Nos1", "Pdyn", "Penk","Npy",
           "Esr1",  
            "Pax6", "Tmem132c", "Dach2", "Lhx2", 
            "Hunk", "Etv1"
            )

pdf(here(plot_dir, "New_DotPlot_marker_genes_Reid_et_al.pdf"), width=6, height=6)
DotPlot_scCustom(seurat, features = genes,  flip_axes=TRUE, colors_use = c("blue", "white", "red")) + NoLegend()
dev.off()

# average expression per cluster
avg <- AverageExpression(seurat, features = genes, return.seurat = FALSE)
avg

min_max_scale <- function(x) (x - min(x)) / (max(x) - min(x))
scaled_mat <- t(apply(avg$SCT, 1, min_max_scale))  # Scale each row

png(here(plot_dir, "Average_expression_per_cluster.png"), width=5, height=4.5, units = "in", res = 300)
ComplexHeatmap::Heatmap(scaled_mat, name = "Expression", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        show_column_names = TRUE, 
                        column_title = "Marker gene expression",
                        #row_title = "Top 5 Marker Genes",
                        # color white to red using colorbewer
                        col = colorRamp2::colorRamp2(c(0, 1), c("white", "red"))
)
dev.off()



# Step 2: Aggregate expression per cluster
# Calculate average expression per cluster for the marker genes
avg_expr <- AverageExpression(seurat, features = genes, return.seurat = FALSE)
expr_matrix <- avg_expr$RNA  # Extract matrix from list

# Step 3: Scale the data for visualization
scaled_expr <- t(scale(t(expr_matrix)))  # Row-wise scaling

# Step 4: Define cluster annotation
cluster_labels <- factor(colnames(scaled_expr))  # Cluster names
cluster_colors <- rainbow(length(unique(cluster_labels)))  # Assign colors
names(cluster_colors) <- unique(cluster_labels)

# Step 5: Generate the heatmap

Heatmap(scaled_expr,
        name = "Expression",
        cluster_rows = FALSE,  # Cluster genes
        cluster_columns = FALSE,  # Don't cluster cell type clusters
        show_column_names = TRUE,
        column_title = "Cell Type Clusters",
        row_title = "Top 5 Marker Genes",
        #col = colorRamp2::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),  # Color scale
        top_annotation = HeatmapAnnotation(Cluster = cluster_labels, 
                                           col = list(Cluster = cluster_colors))
)


# new column order: 15,16, 



