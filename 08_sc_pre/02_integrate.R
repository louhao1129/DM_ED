library(Seurat)
library(qs2)
library(future)
library(string)
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存
rm(list = ls())
setwd("./08_sc_pre")

seurat_filtered = qs_read("./seurat_filtered.qs2")

seurat_filtered <- NormalizeData(seurat_filtered)
seurat_filtered <- FindVariableFeatures(seurat_filtered)
seurat_filtered <- ScaleData(seurat_filtered) # , features = VariableFeatures(seurat_filtered)
seurat_filtered <- RunPCA(seurat_filtered)
ElbowPlot(seurat_filtered, ndims = 50)

# 在metadata中添加分组信息
seurat_filtered@meta.data$group <- ifelse(
  grepl("Diabetes", seurat_filtered@meta.data$orig.ident), "DM",
  ifelse(
    grepl("non-DM", seurat_filtered@meta.data$orig.ident), "non-DM",
    ifelse(
      grepl("Normal", seurat_filtered@meta.data$orig.ident), "Normal", "Other"
    )
  )
)

# 验证分组结果
table(seurat_filtered@meta.data$group)

qs_save(seurat_filtered, "seurat_pca.qs2")

# 整合
# Seurat v5 通过 IntegrateLayers 函数实现简化的整合分析。
# 当前方法支持五种整合方法。这些方法均在低维空间中进行整合，并返回一个降维（即 integrated.{method} ），旨在跨批次共嵌入共享细胞类型

seurat_filtered <- IntegrateLayers(
  object = seurat_filtered, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

seurat_filtered <- FindNeighbors(seurat_filtered, reduction = "harmony", dims = 1:20)
seurat_filtered = FindClusters(
  seurat_filtered,
  resolution = c(seq(0.1, 0.6, 0.1),1)
)

seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:20, reduction = "harmony", reduction.name = "umap")
seurat_filtered <- RunTSNE(seurat_filtered, dims = 1:20, reduction = "harmony", reduction.name = "tsne")

qs_save(seurat_filtered, "seurat_harmony.qs2")

DimPlot(
  seurat_filtered, 
  reduction = "tsne",
  group.by = "RNA_snn_res.0.5",
  split.by = "group"
)

######## 

clustree(seurat_filtered@meta.data, prefix = "RNA_snn_res.")  
ggsave("clusterTree_harmony.png",dpi = 300, width = 10, height = 6)

# 可视化不同分辨率下的聚类UMAP图
for (i in c(seq(0.1, 0.6, 0.1),1)) {
  print(DimPlot(seurat_filtered, reduction = "umap", group.by = paste0("RNA_snn_res.", i)) + labs(title = paste0("resolution: ", i)))
  ggsave(str_glue("umap_harmony_cluster_resolution{i}.png"),dpi=300)
}

# 展示批次效应情况
DimPlot(seurat_filtered, reduction = "umap.harmony", group.by = "group")
ggsave("./02_sc_preprocess/umap_harmony_group.png", dpi=300)