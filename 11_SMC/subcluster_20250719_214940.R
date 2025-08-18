library(Seurat)
library(qs2)
library(future)
library(clustree)
library(patchwork)
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存

SMC = qs_read("SMC.qs2")

SMC <- FindNeighbors(SMC, dims = 1:30, reduction = "pca")
SMC = FindClusters(
  SMC,
  resolution = c(seq(0.1, 1, 0.1)) # 如果设置多个分辨率，那么cluster name也要提供多个或者不提供 
)
# 聚类树
clustree(SMC@meta.data, prefix = "RNA_snn_res.")  
ggsave("clusterTree.png",dpi = 300, width = 10, height = 6)
# TSNE

SMC <- RunTSNE(SMC, dims = 1:20, reduction = "pca", reduction.name = "tsne")
qs_save(SMC, "SMC_cluster.qs2")

#############
p1 = DimPlot(SMC, group.by = "RNA_snn_res.0.3", label = T, label.size = 5, reduction = "tsne")
print(p1)

p2 = DimPlot(SMC, group.by = "RNA_snn_res.0.3", label = T, label.size = 5, split.by = "group", reduction = "tsne")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
print(p2)

p3 = p1 / p2 + plot_layout(guides = "collect")
print(p3)

ggsave("SMC_subcluster.png", dpi=300, height=10, width=10)

#####
SMC = qs_read("./SMC_cluster.qs2")
# 在 excel中书写手动注释结果，然后写入 seurat metadata
Idents(SMC) <- SMC$RNA_snn_res.0.3
annotation_curated_main <- read_excel("./SMC_annotation.xlsx")
new_ids_main <- annotation_curated_main$SMC_type
names(new_ids_main) <- levels(SMC)
SMC <- RenameIdents(SMC, new_ids_main)
SMC@meta.data$SMC_type <- Idents(SMC)

SMC@meta.data |>select(RNA_snn_res.0.3, SMC_type) |> distinct(RNA_snn_res.0.3, .keep_all = TRUE)
DimPlot(SMC, reduction="tsne",group.by = "SMC_type", label = T, label.size = 5, split.by = "group")

DimPlot(SMC, group.by = "SMC_type", reduction="tsne",label = TRUE, label.size = 5) + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
ggsave("SMC_type.png", width = 30, height = 15, dpi=300)

qs_save(SMC, "./SMC_cluster.qs2")
