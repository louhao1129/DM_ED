library(tidyverse)
library(Seurat)
library(qs2)
library(future)
library(viridis)
library(readxl)
library(clustree)
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存
setwd("./09_annotation/")
seurat_harmony_singlet = qs_read("../08_sc_pre/seurat_harmony_singlet.qs2")

seurat_harmony_singlet <- FindNeighbors(seurat_harmony_singlet, reduction = "harmony", dims = 1:20)
seurat_harmony_singlet = FindClusters(
  seurat_harmony_singlet,
  resolution = c(seq(0.1, 1, 0.1))
)

clustree(seurat_harmony_singlet@meta.data, prefix = "RNA_snn_res.")  
ggsave("clusterTree_harmony.png",dpi = 300, width = 10, height = 6)


# Find Marker
Idents(seurat_harmony_singlet) = "RNA_snn_res.0.5"
levels(seurat_harmony_singlet)
cluster.markers <- FindMarkers(seurat_harmony_singlet, ident.1 = 7, ident2=c(0,1,2,6,9,12,13,11,5,17,))
cluster.markers |> 
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(2) & pct.1 >0.8)
write.csv(cluster01.markers, "./05_fibro_annotation/cluster01_AllDiff,csv")

# find markers for every cluster compared to all remaining cells, report only the positive
seurat_harmony_singlet.markers <- FindAllMarkers(seurat_harmony_singlet, only.pos = TRUE)
seurat_harmony_singlet.markers = seurat_harmony_singlet.markers |> 
  filter(pct.1 >0.5)



## 定义marker基因列表
mainmarkers = list(
  "EC" = c("PECAM1", "VWF", "CDH5"),
  "FB" = c("LUM", "COL1A1", "PDGFRA", "DCN"),
  "PC" = c("RGS5", "KCNJ8", "MYOCD","CSPG4"),
  "SMC" = c("ACTA2", "MYH11"),
  "SWC" = c("S100B", "MPZ"),
  "MAC" = c("CD163", "CD68"),
  "T_cells" = c("CD3D", "CD3E", "CD3G")
)

# 查看marker基因的分布，是否特异？一致性如何？
features_list =  c("PECAM1", "VWF", "PTPRC")
FeaturePlot(seurat_harmony_singlet |> subset(RNA_snn_res.0.5 =="15"), features = features_list,reduction = "tsne" ,coord.fixed = T, order = T, cols = viridis(10))
for (i in seq_along(mainmarkers)) {
  p = FeaturePlot(seurat_harmony_singlet, features = mainmarkers[[i]],reduction = "tsne" ,coord.fixed = T, order = T, cols = viridis(10))+
    ggtitle(names(mainmarkers)[i])
  print(p)
  ggsave(paste0("FeaturePlot_mainmarkers_", names(mainmarkers)[i], ".png"), path = "feature_plot", width = 10, height = 10, units = "cm")
}

# 尝试多种分辨率
# 结合聚类树查看分群演化情况
# 对于混合的类型，先标注xx/xx混合，可以提取子集再次标注化、聚类
DotPlot(seurat_harmony_singlet, features = mainmarkers, group.by = "RNA_snn_res.0.5") + 
  # coord_flip() + 
  scale_color_viridis()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  ggtitle("res = 0.5")
ggsave("DotPlot_mainmarkers.png", width = 15, height = 10, bg = "white")

DimPlot(seurat_harmony_singlet, reduction = "tsne", group.by = "RNA_snn_res.0.5")
FeaturePlot(seurat_harmony_singlet, features = mainmarkers[[1]], coord.fixed = T, order = T, cols = viridis(10))

############ 将注释结果写入metadata
# 在 excel中书写手动注释结果，然后写入 seurat metadata
Idents(seurat_harmony_singlet) <- seurat_harmony_singlet$RNA_snn_res.0.5
annotation_curated_main <- read_excel("./main_annotation.xlsx")
new_ids_main <- annotation_curated_main$main_cell_type
names(new_ids_main) <- levels(seurat_harmony_singlet)
seurat_harmony_singlet <- RenameIdents(seurat_harmony_singlet, new_ids_main)
seurat_harmony_singlet@meta.data$main_cell_type <- Idents(seurat_harmony_singlet)

# 查看注释情况
seurat_harmony_singlet@meta.data |>select(RNA_snn_res.0.5, main_cell_type) |> distinct(RNA_snn_res.0.5, .keep_all = TRUE)
DimPlot(seurat_harmony_singlet, reduction="tsne",group.by = "main_cell_type", label = T, label.size = 5)

DimPlot(seurat_harmony_singlet, group.by = "main_cell_type", reduction="tsne",label = TRUE, label.size = 5) + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
ggsave("main_cell_type.png", width = 30, height = 15, dpi=300)

qs_save(seurat_harmony_singlet, "./seurat_annotation.qs2")
