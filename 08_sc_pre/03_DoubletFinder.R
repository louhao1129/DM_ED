library(Seurat)
library(qs2)
library(future)
library(DoubletFinder) # https://github.com/chris-mcginnis-ucsf/DoubletFinder?tab=readme-ov-file
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存
rm(list = ls())
setwd("./08_sc_pre")

seurat_harmony = qs_read("./seurat_harmony.qs2")
seurat_harmony <- JoinLayers(seurat_harmony)

# DoubletFinder有几个比较重要的参数和概念：
# 1）pN：指生成的模拟双细胞的数量，以合并的真实-模拟数据的比例表示，一般这个值设置为25%
# 2）pK：指最近邻细胞之间比例。没有设置默认值，pK需要根据不同的数据进行调整，寻找最优pK也是DoubletFinder中很重要的一步。
# 3）Homotypic：同源双细胞的比例，是指由两个或多个相同类型的细胞融合而成的假阳性细胞，而我们真正需要去除的是异源双细胞。
# 4）nExp：表示预计的异源双细胞数量

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
# # 对"scRNA_harmony"这个单细胞对象进行pN-pK参数扫描，以生成人工双细胞并计算每个细胞的pANN值
sweep.res.list <- paramSweep(seurat_harmony, PCs = 1:20, sct = FALSE)
# 对参数扫描的结果进行汇总，计算每种pN和pK组合下的双细胞检测指标。这些指标可以用来评估不同参数下的双细胞检测效果，并选择最优的参数。参数GT表示是否提供了真实的双细胞标签，此处没有提供
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
# 根据汇总结果找到最优的pK参数。
bcmvn<- find.pK(sweep.stats)
# 提取出全局最优的pK值，储存于"pK_bcmvn",最高点所对应的pK就是全局最优pK
pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
# 估计的同源双细胞（由两个或多个相同类型的细胞融合而成的假阳性细胞，它们通常比异源双细胞更难以检测和去除）的比例     
homotypic.prop <- modelHomotypic(seurat_harmony@meta.data$RNA_snn_res.0.5)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
# 计算总的双细胞数量（假设双细胞形成率为 4.6%）
# 我预期的双倍费率是多少？答案：这取决于您的平台（10x、parse 等）
nExp_poi <- round(0.046*nrow(seurat_harmony@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
# 使用确定好的参数鉴定doublets
# DoubletFinder 性能在很大程度上是 pN 不变的,默认值设置为 25%
plan("multisession", workers =2)
options(future.globals.maxSize = 60*1024^3) # 每个线程分配运行内存
seurat_harmony <- doubletFinder(seurat_harmony, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
seurat_harmony <- doubletFinder(seurat_harmony, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_17_2990", sct = FALSE)

# 可视化
DimPlot(seurat_harmony, reduction = "tsne", group.by = "DF.classifications_0.25_17_2513", raster = FALSE)
DimPlot(seurat_harmony, reduction = "tsne", group.by = "RNA_snn_res.0.5")
qs_save(seurat_harmony,"./seurat_harmony_doublefinder.qs2")
# 将双细胞剔除后生成新的对象，便于我们后续分析
seurat_harmony_singlet <- subset(seurat_harmony_doublefinder, subset = DF.classifications_0.25_17_2513== "Singlet")
DimPlot(seurat_harmony_singlet, reduction = "tsne", group.by = "RNA_snn_res.0.5", split.by = "group")

qs_save(seurat_harmony_singlet, "./seurat_harmony_singlet.qs2")
rm(list = ls())
