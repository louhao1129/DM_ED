rm(list = ls())
library(Seurat)
library(qs2)
library(future)
library(readr)
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(UCell)
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存
setwd("./10_signature/")

hub_gene = read_csv("../07_cytoscape/net.network.txt_Degree_top15node.csv")

seurat_annotation = qs_read("../09_annotation/seurat_annotation.qs2")

all_list = hub_gene |> pull(name) |> as.list()
down_list = hub_gene |> filter(Type == "down") |> pull(name) |> as.list()
up_list = hub_gene |> filter(Type == "up") |> pull(name) |> as.list()

marker_score_all <- AddModuleScore_UCell(seurat_annotation,
                                     features=all_list)
marker_score_down <- AddModuleScore_UCell(seurat_annotation,
                                     features=down_list)

marker_score_up <- AddModuleScore_UCell(seurat_annotation,
                                     features=up_list)

p_down = FeaturePlot(marker_score_down , features =  "signature_3_UCell",reduction = "tsne" ,coord.fixed = T, order = T, split.by = "group")
p_up = FeaturePlot(marker_score_up , features =  "signature_12_UCell",reduction = "tsne" ,coord.fixed = T, order = T, split.by = "group")
p_all = FeaturePlot(marker_score_all , features =  "signature_15_UCell",reduction = "tsne" ,coord.fixed = T, order = T, split.by = "group")

p_ucell = p_up / p_down / p_all
p_ucell + plot_annotation(
  title = "signature ucell score"
)
ggsave("signature_ucell_score.png", p_ucell, dpi=300)


# seurat_annotation = AddModuleScore(
#   seurat_annotation,
#   up_list,
#   name = "up"
# )

# seurat_annotation = AddModuleScore(
#   seurat_annotation,
#   down_list,
#   name = "down"
# )



# p1 = FeaturePlot(seurat_annotation , features = "up12",reduction = "tsne" ,coord.fixed = T, order = T, split.by = "group")
# p2 = FeaturePlot(seurat_annotation , features = "down3",reduction = "tsne" ,coord.fixed = T, order = T, split.by = "group")
# p1
# p2

# # for (i in seq_along(down_list)) {
# #   p = FeaturePlot(seurat_annotation, features = down_list[[i]],reduction = "tsne" ,coord.fixed = T, order = T, cols = viridis(10))+
# #     ggtitle(names(down_list)[i])
# #   print(p)
# #   #ggsave(paste0("FeaturePlot_mainmarkers_", names(mainmarkers)[i], ".png"), path = "feature_plot", width = 10, height = 10, units = "cm")
# # }

# p3 = DimPlot(seurat_annotation, reduction="tsne",group.by = "main_cell_type", label = T, label.size = 5, split.by = "group")
# p3

# p4 = p3 / p2 / p1 + plot_layout(guides = "collect")

# ggsave("./signature.png", p4, dpi=300, bg="white")


# DimPlot(seurat_annotation, reduction="tsne",group.by = "main_cell_type", label = T, label.size = 5)
# ggsave("./main_cell_type_tsne.png", dpi=300)
