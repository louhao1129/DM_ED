# if (!requireNamespace("irGSEA", quietly = TRUE)) { 
#     devtools::install_github("chuiqin/irGSEA")
# }
library(tidyverse)
library(Seurat)
library(qs2)
library(future)
library(viridis)
library(UCell)
library(GSEABase)
library(irGSEA) # https://www.jianshu.com/p/463dd6e2986f
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存
rm(list = ls())
setwd("./13_irGSEA/")
SMC = qs_read("../11_SMC/SMC_cluster.qs2")
Idents(SMC) = "RNA_snn_res.0.1"
levels(SMC)

# geneset=read.gmt("./KEGG_metabolism_nc.gmt") # 其实就是一个dataframe
# length(unique(geneset$term)) # 查看通路数量
# table(geneset$term) # 查看每个通路的基因数量

gmt_file <- "./KEGG_metabolism_nc.gmt"
gene_sets <- getGmt(gmt_file)
# 转换为 list（irGSEA 接受 list）
custom_gene_list <- geneIds(gene_sets)

SMC.final <- irGSEA.score(
  object = SMC, 
  assay = "RNA", 
  slot = "data", 
  seeds = 123, 
  ncores = 10,
  min.cells = 3, 
  min.feature = 0,
  custom = TRUE,
  geneset = custom_gene_list,   # ← 用你的基因集
  msigdb = FALSE,               # ← 关闭默认 msigdb
  species = "Homo sapiens",
  method = c("AUCell", "UCell", "singscore", "ssgsea"),
  kcdf = 'Gaussian'
)

qs_save(SMC.final, "./SMC_irGSEA.qs2")

# 返回一个Seurat对象，富集分数矩阵存放在RNA外的assay中
Seurat::Assays(SMC.final)
#> [1] "RNA"       "AUCell"    "UCell"     "singscore" "ssgsea"
result.dge <- irGSEA.integrate(object = SMC.final, 
                               group.by = "RNA_snn_res.0.1",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))

qs_save(result.dge, "deg_irGSEA.qs2")
###### 可视化
# 全局展示
# 热图展示了综合评价中具体基因集在每个细胞亚群是否具有统计学意义差异；
# 其中，浅蓝色的格子无统计学差异，红色的格子具有统计学差异。格子中的星号越多，格子的P值越小；
# 左边的聚类树代表不同基因集在不同细胞亚群中表达模式的相似性；
# 上方的条形图分别代表不同的细胞亚群，以及差异基因集在细胞亚群中是呈现上调还是下调趋势；
# 你还可以把method从'RRA"换成“ssgsea”，展示特定基因集富集分析方法中差异上调或差异下调的基因集；
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.heatmap.plot
ggsave("irGSEA.heatmap.plot.png", dpi=300, plot=irGSEA.heatmap.plot, height=10, width=20)

# 气泡图展示了综合评价中具体基因集在每个细胞亚群是否具有统计学意义差异；
# 其中，浅蓝色的点无统计学差异，红色的点具有统计学差异。点越大，P值越小；
# 左边的聚类树代表不同基因集在不同细胞亚群中表达模式的相似性；
# 上方的条形图分别代表不同的细胞亚群，以及差异基因集在目标细胞亚群中是呈现上调还是下调趋势；
# 你还可以把method从'RRA"换成“ssgsea”，展示特定基因集富集分析方法中差异上调或差异下调的基因集；
irGSEA.bubble.plot <- irGSEA.bubble(object = result.dge, 
                                    method = "RRA", 
                                    top = 50)
irGSEA.bubble.plot
ggsave("irGSEA.bubble.plot.png", dpi=300, plot=irGSEA.bubble.plot)

# upset图展示了综合评估中每个细胞亚群具有统计学意义差异的基因集的数目，以及不同细胞亚群之间具有交集的差异基因集数目；
# 左边不同颜色的条形图代表不同的细胞亚群；
# 上方的条形图代表具有交集的差异基因集的数目；
# 中间的气泡图单个点代表单个细胞亚群，多个点连线代表多个细胞亚群取交集.；
irGSEA.upset.plot <- irGSEA.upset(object = result.dge, 
                                  method = "RRA")
irGSEA.upset.plot

# 堆叠柱状图具体展示每种基因集富集分析方法中每种细胞亚群中上调、下调和没有统计学差异的基因集数目；
# 上方的条形代表每个亚群中不同方法中差异的基因数目，红色代表上调的差异基因集，蓝色代表下调的差异基因集；
# 中间的柱形图代表每个亚群中不同方法中上调、下调和没有统计学意义的基因集的比例；
irGSEA.barplot.plot <- irGSEA.barplot(object = result.dge,
                                      method = c("AUCell", "UCell", "singscore",
                                                 "ssgsea"))
irGSEA.barplot.plot


# 局部展示

# 密度散点图将基因集的富集分数和细胞亚群在低维空间的投影结合起来，
# 展示了特定基因集在空间上的表达水平。其中，颜色越黄，代表富集分数越高；
names(custom_gene_list)
scatterplot <- irGSEA.density.scatterplot(object = SMC.final,
                             method = "UCell",
                             show.geneset = "Glycolysis / Gluconeogenesis",
                             reduction = "tsne")
scatterplot
ggsave("irGSEA.scatterplot.plot.png", dpi=300, plot=scatterplot)

# 半小提琴图同时以小提琴图（左边）和箱线图（右边）进行展示。不同颜色代表不同的细胞亚群；
halfvlnplot <- irGSEA.halfvlnplot(object = SMC.final,
                                  method = "UCell",
                                  show.geneset = "Glycolysis / Gluconeogenesis")
halfvlnplot

# 山峦图中上方的核密度曲线展示了数据的主要分布，下方的条形编码图展示了细胞亚群具体的数量。
# 不同颜色代表不同的细胞亚群，而横坐标代表不同的表达水平；
ridgeplot <- irGSEA.ridgeplot(object = SMC.final,
                              method = "UCell",
                              show.geneset = "Glycolysis / Gluconeogenesis")
ridgeplot
