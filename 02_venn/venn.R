library(VennDiagram)      
library(dplyr)
setwd("./02_venn/")
rm(list = ls())

listFile="interGenes.List.txt"        # 输出交集基因的列表文件
typeFile="interGenes.Type.txt"        # 输出交集基因的属性文件
upList=list()        # 上调基因的列表
downList=list()      # 下调基因的列表

# 获得上调基因和下调基因的列表，用于绘制韦恩图
deg1 = read.csv("../01_bulk_anno_diff/GSE10804_deg.csv")
deg2 = read.csv("../01_bulk_anno_diff/GSE25724_deg.csv")

deg1_up = deg1 |> filter(abs(logFC) >0.5 & P.Value <0.05 & logFC>0)|> pull(symbol)
deg1_down = deg1 |> filter(abs(logFC) >0.5 & P.Value <0.05 & logFC<0) |> pull(symbol)
deg2_up = deg2 |> filter(abs(logFC) >0.5 & P.Value <0.05 & logFC>0) |> pull(symbol)
deg2_down = deg2 |> filter(abs(logFC) >0.5 & P.Value <0.05 & logFC<0) |> pull(symbol)

upList = list(deg1_up, deg2_up)
names(upList) = c("GSE10804_up", "GSE25724_up")
downList = list(deg1_down, deg2_down)
names(downList) = c("GSE10804_down", "GSE25724_down")

# 绘制上调基因韦恩图
venn.plot=venn.diagram(upList,filename=NULL,fill=c("cornflowerblue", "darkorchid1"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.1)
pdf(file="up.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#绘制下调基因韦恩图
venn.plot=venn.diagram(downList,filename=NULL,fill=c("cornflowerblue", "darkorchid1"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.1)
pdf(file="down.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

# 输出交集基因文件
upInterGenes=Reduce(intersect, upList)        #交集上调基因,Reduce(intersect, list(v1, v2, v3))
upTab=cbind(upInterGenes, "up")
downInterGenes=Reduce(intersect, downList)    #交集下调基因
downTab=cbind(downInterGenes, "down")
interTab=rbind(upTab, downTab)
colnames(interTab)=c("Gene", "Type") # 合并
#输出交集基因的属性文件
write.table(file=typeFile, interTab, sep="\t", quote=F, col.names=T, row.names=F)
#输出交集基因的列表文件
write.table(file=listFile, interTab[,1], sep="\t", quote=F, col.names=F, row.names=F)

