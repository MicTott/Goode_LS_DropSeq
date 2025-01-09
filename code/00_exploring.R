library(SingleCellExperiment)
library(Seurat)
library(here)
library(scran)
library(scater)

load(here("processed-data", "DLSclusters.Robj"))
seurat.obj

seurat.v3 <- UpdateSeuratObject(seurat.obj)
sce <- as.SingleCellExperiment(seurat.v3)
sce

colnames(colData(sce))
# [1] "orig.ident"   "res.1"        "nCount_RNA"   "nFeature_RNA" "ident"  

unique(sce$orig.ident)

plotReducedDim(sce,dimred="TSNE", colour_by="ident")
plotReducedDim(sce,dimred="TSNE", colour_by="orig.ident")



# get mito genes
is.mito <- grepl("^mt-", rownames(sce))

# subset to sce to is mito
sce.mito <- sce[rownames(sce) %in% is.mito,]

# add scuttle qc metrics
sce <- scuttle::addPerCellQC(sce, subsets=list(mito=is.mito))

plotColData(sce,x="orig.ident", y="total", colour_by="orig.ident")
plotColData(sce,x="orig.ident", y="detected", colour_by="orig.ident")
plotColData(sce,x="orig.ident", y="subsets_mito_percent", colour_by="orig.ident")


# ====== Add discards by MADs ========
sce$sum_discard <- isOutlier(sce$sum, nmads=3, type="higher")
sce$detected_discard <- isOutlier(sce$detected, nmads=3, type="higher")
sce$subsets_mito_percent_discard <- isOutlier(sce$subsets_mito_percent, nmads=3, type="higher")


plotColData(sce,x="orig.ident", y="sum", colour_by="sum_discard")
plotColData(sce,x="orig.ident", y="detected", colour_by="detected_discard")
plotColData(sce,x="orig.ident", y="subsets_mito_percent", colour_by="subsets_mito_percent_discard")

sce$sample_id <- sce$orig.ident
# drop discards
sce <- sce[,!sce$sum_discard]
sce <- sce[,!sce$detected_discard]
sce <- sce[,!sce$subsets_mito_percent_discard]



# ===== Normalization =====
set.seed(100)
clust <- quickCluster(sce, BPPARAM=BiocParallel::MulticoreParam(8)) 
sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1)
sce <- logNormCounts(sce)



# ===== Feature selection =====
dec <- modelGeneVar(sce)
chosen <- getTopHVGs(dec, n=500)


# ===== Dimensionality Reduction =====
set.seed(255) # See below.
sce <- runPCA(sce, subset_row=chosen) 

percent.var <- attr(reducedDim(sce,"PCA"), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")


set.seed(376)
sce <- runUMAP(sce, dimred="PCA", pca=20, min_dist=0.3)
plotReducedDim(sce, dimred="UMAP", colour_by="sample_id")


# ===== Clustering =====
library(bluster)

set.seed(492)
clust.leiden <- clusterCells(sce, use.dimred="PCA", 
                              BLUSPARAM=NNGraphParam(cluster.fun="louvain"))

sce$louvain <- clust.leiden

plotReducedDim(sce, dimred="UMAP", colour_by="louvain", text_by="louvain")
plotReducedDim(sce, dimred="UMAP", colour_by="Slc17a6", text_by="louvain")


# ===== Marker gene plotting =====
features <- c("Mbp", "Mobp",
              "Snap25", "Syt1",
              "Slc17a7", "Slc17a6",
              "Slc32a1", "Slc6a1",
              "Gad1", "Gad2"
              )

plotDots(sce, features, group="louvain")



# drop 1, 2, 13, 14, 16, 24, 26
sce <- sce[,!sce$louvain %in% c(1, 2, 13, 14, 16, 24, 26)]

# ===== Save as rds =====
saveRDS(sce, here("processed-data", "sce_only_GABA.rds"))



# ========= Sub-clustering GABAergic cells ==========

# ===== Normalization =====
clust <- quickCluster(sce, BPPARAM=BiocParallel::MulticoreParam(8)) 
sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1)
sce <- logNormCounts(sce)



# ===== Feature selection =====
dec <- modelGeneVar(sce)
chosen <- getTopHVGs(dec, n=15000)


# ===== Dimensionality Reduction =====
sce <- runPCA(sce, subset_row=chosen) 

percent.var <- attr(reducedDim(sce,"PCA"), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")


sce <- runUMAP(sce, dimred="PCA", pca=50, min_dist=0.3)
plotReducedDim(sce, dimred="UMAP", colour_by="sample_id")
plotReducedDim(sce, dimred="UMAP", colour_by="Npy")
