library(limma) # 芯片数据差异分析，使用limma包
library(GEOquery)
library(dplyr)
library(tidyverse)
rm(list = ls())

load("./GSE95849_eSet.Rdata")
gse = gset[[1]]
pdata = pData(gse)
group_list = pdata$title
group_list = group_list[1:12]


exp.4 = readRDS("./GSE95849_exp.rds")
exp.4 = exp.4[,1:12]
boxplot(exp.4)
##除去批次效应，同一个实验做多次重复，结果不大可能一致
exp.5=normalizeBetweenArrays(exp.4) # 同一芯片（平台）批次效应去除
boxplot(exp.5)
# 进行normalize之后，比较才有可比性，才能进行差异分析
a = 1:6
a2 = 7:12
exp.5=exp.5[,c(a2,a)] # 改变排列顺序，normal在前
saveRDS(exp.5, "exp.5.rds")

##这里我们创建分组信息
group_list=c(rep('control',6),rep('T',6))
##转化为因子,为了后面差异分析的矩阵构建
group_list=factor(group_list)
## 强制限定顺序
group_list <- relevel(group_list, ref="control") # ref=control，代表control为参考组，在前，和表达矩阵中样本顺序相统一

##构建差异分析的矩阵
design=model.matrix(~ group_list)
# 分组矩阵有两中构建方法，我们这里~前面没有0，说明第一列作为比对列（control），第二列与第一列相比
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp.5)


##lmFit()：线性拟合模型构建
fit=lmFit(exp.5,design)
##eBayes()使用trend=TRUE对标准误差进行经验贝叶斯平滑，计算每个对比中每个基因的moderated t-statistic和log-odds。
fit=eBayes(fit) 
##topTable()给出一个最有可能在给定对比下差异表达的基因列表。
allDiff=topTable(fit,coef=2,adjust='fdr',number=Inf) # coef=2代表分组矩阵第二列与第一列比较的结果
##topTable函数的coef参数，coef=2是指design的第2列，即tumour，即把tumour与normal进行对比
allDiff$symbol = rownames(allDiff)
write.table(allDiff,file = "allDiff.txt",sep = "\t",col.names = NA)
write.csv(allDiff, "allDIff.csv")
saveRDS(allDiff, "GSE95849_deg.rds")
