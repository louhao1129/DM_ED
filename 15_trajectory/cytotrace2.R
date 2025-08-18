library(CytoTRACE2)
library(Seurat)
library(tidyverse)
library(patchwork)
library(patchwork)
library(qs2)

setwd("./15_trajectory/")
SMC = qs_read("../11_SMC/SMC_cluster.qs2")
Idents(SMC) = "RNA_snn_res.0.1"
levels(SMC)
SMC <- RunUMAP(SMC, dims = 1:20)
DimPlot(SMC, reduction = "umap", group.by = "RNA_snn_res.0.1")

cytotrace2_result = cytotrace2(
  SMC,
  species = "human",
  is_seurat = TRUE, # 接受一个Seurat对象，返回一个Seurat对象
  slot_type = "counts",
  batch_size = 10000,
  smooth_batch_size = 1000,
  parallelize_models = TRUE,
  parallelize_smoothing = TRUE,
  ncores = NULL, # Windows OS can run only on 1 core
  seed = 14
)
qs_save(cytotrace2_result, "SMC_ctyotrace.qs2")
# absolute developmental potential of each cell, which we term as "potency score", 
# as a continuous value ranging from 0 (differentiated) to 1 (stem cells capable of generating an entire multicellular organism). 
annotation <- data.frame(phenotype = SMC@meta.data$RNA_snn_res.0.1) %>% 
  set_rownames(., colnames(SMC)) # 把 annotation 这个 data.frame 的行名设置为 SMC 对象中细胞的名字（即 colnames(SMC)）

plots <- plotData(
  cytotrace2_result,
  annotation = annotation,
  expression_data = NULL,
  is_seurat = TRUE,
  slot_type = "counts",
  pc_dims = 30,
  seed = 14
)

plots$CytoTRACE2_UMAP
plots$CytoTRACE2_Relative_UMAP
plots$CytoTRACE2_Boxplot_byPheno

FeaturePlot(cytotrace2_result, reduction = "umap", features = "CytoTRACE2_Score")
DimPlot(cytotrace2_result, reduction = "tsne", group.by = "RNA_snn_res.0.1")
DimPlot(cytotrace2_result, reduction = "umap", group.by = "RNA_snn_res.0.1")


#### 统计检验
df = cytotrace2_result@meta.data
# 分组 0 vs 1
group_0_vs_1 <- df[df$RNA_snn_res.0.1 %in% c(0, 1), ]
test_0_vs_1 <- wilcox.test(CytoTRACE2_Score ~ RNA_snn_res.0.1, data = group_0_vs_1)

# 分组 0 vs 2
group_0_vs_2 <- df[df$RNA_snn_res.0.1 %in% c(0, 2), ]
test_0_vs_2 <- wilcox.test(CytoTRACE2_Score ~ RNA_snn_res.0.1, data = group_0_vs_2)

# 输出结果
test_0_vs_1
# W = 1414468, p-value < 2.2e-16
test_0_vs_2
# W = 528481, p-value = 0.06817
