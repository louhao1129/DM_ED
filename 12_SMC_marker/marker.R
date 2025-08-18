library(tidyverse)
library(Seurat)
library(qs2)
library(future)
library(viridis)
library(readxl)
library(clustree)
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存

SMC = qs_read("../11_SMC/SMC_cluster.qs2")
# Find Marker
Idents(SMC) = "RNA_snn_res.0.1"
levels(SMC)
cluster.markers <- FindMarkers(SMC, ident.1 = 0)
cluster.markers |> 
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5) & pct.1 >0.7 & pct.2<0.3)
write.csv(cluster01.markers, ".")

# find markers for every cluster compared to all remaining cells, report only the positive
SMC.markers <- FindAllMarkers(SMC, only.pos = TRUE)
SMC.markers |> 
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5 & pct.1 >0.6 & pct.2<0.4)
write.csv(SMC.markers, "SMC_cluster_marker.csv")

## 定义marker基因列表

SMC_marker = list(
  "cluster0"=c("S100A4", "IFITM2", "PLAC9"),
  "cluster1"=c("MTRNR2L12","AQP1", "HSPB8", "PPP1R1A"),
  "cluster2"=c("CES1", "CTGF", "KCNAB1", "PCP4")
)
# 不要被贩卖焦虑啦，放轻松点呀，嘻嘻，要的本来就不多嘛
# 查看marker基因的分布
features_list =  c("PECAM1", "VWF", "PTPRC")
FeaturePlot(SMC |> subset(RNA_snn_res.0.5 =="15"), features = features_list,reduction = "tsne" ,coord.fixed = T, order = T, cols = viridis(10))
for (i in seq_along(SMC_marker)) {
  p = FeaturePlot(SMC, features = SMC_marker[[i]],reduction = "tsne" ,coord.fixed = T, order = T, cols = viridis(10))+
    ggtitle(names(SMC_marker[[i]]))
  print(p)
  ggsave(paste0("FeaturePlot_SMC_marker_", names(SMC_marker)[i], ".png"), width = 5, height = 5, units = "cm")
}


# 气泡图
DotPlot(SMC, features = SMC_marker, group.by = "RNA_snn_res.0.1") + 
  # coord_flip() + 
  scale_color_viridis()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  ggtitle("res = 0.1")
ggsave("DotPlot_SMC_markers.png", width = 8, height = 5, bg = "white")

DimPlot(SMC, reduction = "tsne", group.by = "RNA_snn_res.0.1", split.by = "group")
ggsave("SMC_subcluster.png", dpi=300, width=10, height=5)

FeaturePlot(SMC, features = SMC_marker[[3]], reduction="tsne",coord.fixed = T, order = T, cols = viridis(10))
ggsave("FeaturePlot_marker_cluster2.png", dpi=300)
