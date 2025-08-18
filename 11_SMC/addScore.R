library(viridis)
library(UCell)
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存
rm(list = ls())
setwd("./11_SMC/")

hub_gene = read_csv("../07_cytoscape/net.network.txt_Degree_top15node.csv")
SMC = qs_read("SMC_cluster.qs2")

down_list = hub_gene |> filter(Type == "down") |> pull(name) |> as.list()
up_list = hub_gene |> filter(Type == "up") |> pull(name) |> as.list()

SMC = AddModuleScore(
  SMC,
  down_list,
  name = "down"
)
# 我这可是疾病组下调的基因啊！！！
FeaturePlot(SMC , features = "down3",reduction = "tsne" ,coord.fixed = T, order = T, split.by = "group")
DimPlot(SMC, reduction="tsne",group.by = "RNA_snn_res.0.3", label = T, label.size = 5, split.by = "group")

# 或许应该考虑从单细胞的角度出发做共病，先发现共病细胞群