library(WGCNA)
library(qs2)
library(tidyverse)
library(GEOquery)
options(stringsAsFactors = FALSE)

# https://space.bilibili.com/3546800304687494
getwd()
setwd("./DM_ED/01_bulk_GSE25724/")
# 创建目录（若存在则静默跳过）
fs::dir_create("WGCNA_output", recurse = TRUE)

# 表达矩阵清洗
mrna_expr = readRDS("./GSE25724_exp.rds")
head(mrna_expr[1:4, 1:4]) #是否取log？或者做过标准化
#mrna_expr = log2(mrna_expr+1)

# 挑选基因，采用绝对中位差mad过滤数据
mrna_exp_mad = apply(mrna_expr, 1, mad)
mrna_exp_mad = order(mrna_exp_mad, decreasing = TRUE)
mrna_exp_num = mrna_exp_mad[1:5000]
mrna_exp_filter = mrna_expr[mrna_exp_num,]
datExpr0 = t(mrna_exp_filter)

# 对样本和基因进行评估，剔除低质量的样本和基因
gsg = goodSamplesGenes(datExpr0, verbose = 3) 
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# 聚类树，对离群样本进行剔除
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = paste0("WGCNA_output/", "1_sampleClustering.pdf"), width = 25, height = 12)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 110000, col = "red")# 样本剪切高度，在这个高度正好把离群样本和其他样本分割开来
dev.off()

#选做，根据高度过滤样本
# Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = 110000, minSize = 10)
# print(table(clust)) # 0表示要剔除的样本，1表示要保留的样本

# # clust 1 contains the samples we want to keep.
# keepSamples = (clust==1)
# datExpr0 = datExpr0[keepSamples, ]

qs_save(datExpr0, "WGCNA_output/datExpr0.qs2")
rm(list = ls())


##### 表型信息处理
load("./GSE25724_eSet.Rdata")
datExpr0 = qs_read("WGCNA_output/datExpr0.qs2")

# 数据分组
gse = gset[[1]]
pd = pData(gse)
glimpse(pd)
Group=ifelse(str_detect(pd$title,"Non"), # 自行寻找分组信息所在的组合分组的情况
             "Control",
             "diabetes")
Group = data.frame(
  sample = rownames(datExpr0),
  Group = Group
)
traitdata = Group |>
  mutate(diabetes = ifelse(Group == "diabetes", 1, 0)) |> 
  column_to_rownames("sample") |> 
  select(diabetes)
qs_save(traitdata, "WGCNA_output/traitdata.qs2")
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(traitdata, signed = FALSE)
#sizeGrWindow(12,12)
pdf(file=paste0("WGCNA_output/", "2_Sample dendrogram and trait heatmap.pdf"),width=25,height=10)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitdata),
                    main = "Sample dendrogram and trait heatmap",marAll = c(3,10,3,3))
dev.off()



rm(list = ls())

# 逐步法识别基因模块
traitdata = qs_read("WGCNA_output/traitdata.qs2")
datExpr0 = qs_read("WGCNA_output/datExpr0.qs2")
# 多线程
enableWGCNAThreads(nThreads = 8)
# Choose a set of soft-thresholding powers
powers = c(1:20)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results:
#sizeGrWindow(9, 5)
pdf(file=paste0("WGCNA_output/", "3_Scale independence.pdf"),width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#chose the softPower
softPower =sft$powerEstimate # 提取最佳软阈值
print(sft$powerEstimate)
adjacency = adjacency(datExpr0, power = softPower) # 基因表达矩阵 -> 相关性矩阵 -> 邻接矩阵

k = as.vector(apply(adjacency, 2, sum, na.rm=T))

par(mfrow = c(1,2))
cex = 0.9
hist(k)
scaleFreePlot(k, main = "check scale free topology") # 注意R^2是否大于0.85

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency) # 邻接矩阵 -> 拓扑矩阵
dissTOM = 1-TOM # 计算基因相异度

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)

#sizeGrWindow(12,9)
pdf(file=paste0("WGCNA_output/", "4_Gene clustering on TOM-based dissimilarity.pdf"),width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# 动态剪切树逐步识别网络模块
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, # 基因聚类结果
                            distM = dissTOM, # 基因相异度
                            deepSplit = 2, # 调整划分模块的敏感度，数值（0-4）越大越敏感，得到的模块越多
                            pamRespectsDendro = FALSE, # 不对相似度高的模块进行自动合并
                            minClusterSize = minModuleSize)  # 模块包含的最小基因数量，数值越小，小的模块会被保留，模块越多
table(dynamicMods) # 模块数和每个模块所含的基因数

# 绘制基因聚类树图
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file=paste0("WGCNA_output/", "5_Dynamic Tree Cut.pdf"),width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes,计算每个模块的特征基因表达量
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes # 每个基因模块在每个样本中的表达量
# Calculate dissimilarity of module eigengenes，计算模块特征基因相异度
MEDiss = 1-cor(MEs);
# Cluster module eigengenes，根据相似度进行聚类
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file=paste0("WGCNA_output/", "6_Clustering of module eigengenes.pdf"),width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.5 # 剪切高度可修改，低于高度的模块可以进行合并
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function，合并相异度较低的模块
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors，提取合并后模块对应的颜色
mergedColors = merge$colors
print(length(table(mergedColors))) # 查看模块数量
# Eigengenes of the new merged modules，提取合并之后每个基因模块的表达量
mergedMEs = merge$newMEs

# 可视化合并之后的网络模块
#sizeGrWindow(12, 9)
pdf(file=paste0("WGCNA_output/", "7_merged dynamic.pdf"), width = 9, height = 6.5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
save(softPower, MEs, moduleLabels, moduleColors, geneTree, file = "WGCNA_output/逐步法识别基因模块.rdata")

rm(list = ls())

####### 表型信息关联分析
datExpr0 = qs_read("WGCNA_output/datExpr0.qs2")
datTraits = qs_read("WGCNA_output/traitdata.qs2")
load("WGCNA_output/逐步法识别基因模块.rdata")
# 对模块特征基因表达矩阵进行排序，使相似的基因模块（通过相关性测量）相邻
MEs = orderMEs(MEs)
names = rownames(MEs)
# datTraits = datTraits[names,] # 样本对其
datTraits <- datTraits[names, , drop = FALSE]
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
# 计算模块和表型信息相关性
moduleTraitCor = cor(MEs, datTraits, use = "p") # 行名是模块名，列名是表型信息
# 计算模块和表型信息相关性p值
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(10,6)
pdf(file=paste0("WGCNA_output/", "8_Module-trait relationships.pdf"),width=9,height=8)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,font.lab.x = 0.5,font.lab.y = 0.5,
  zlim = c(-1,1),
  main = paste("Module-trait relationships"))
dev.off()

# 可以把相关性最高的模块提取出来做最后的分析


##### 网络模块可视化
# 提取感兴趣模块中的基因
modules = c("antiquewhite4")
inModule = (moduleColors==modules) # 基因所在位置
modGenes = colnames(datExpr0)[inModule]
write.csv(modGenes, "modGenes.csv", row.names=F)
