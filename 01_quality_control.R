library(Seurat)
library(here)

load(here("processed-data", "DLSclusters.Robj"))
seurat.obj


colnames(seurat@meta.data)
# [1] "nCount_RNA"      "nFeature_RNA"    "percent.mito"    "Sample_ID"      
# [5] "nCount_SCT"      "nFeature_SCT"    "SCT_snn_res.0.5" "seurat_clusters"
# [9] "broad_celltype" 

# violins of Qc metrics
pdf(here("plots", "01_quality_control", "qc_violins.pdf"), width=12, height=4)
VlnPlot(seurat, group.by="Sample_ID", features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3, pt.size=0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()
