library(ggVennDiagram)
library(ggplot2)
setwd("./02_ggvenn/")

listFile="interGenes.List.txt"        # 输出交集基因的列表文件
typeFile="interGenes.Type.txt"        # 输出交集基因的属性文件
upList=list()        # 上调基因的列表
downList=list()      # 下调基因的列表

deg1 = read.csv("../01_bulk_GSE10804/GSE10804_deg.csv")
deg2 = read.csv("../01_bulk_GSE25724/GSE25724_deg.csv")

# 获得上调基因和下调基因的列表，用于绘制韦恩图
deg1 = read.csv("../01_bulk_anno_diff/GSE10804_deg.csv")
deg2 = read.csv("../01_bulk_anno_diff/GSE25724_deg.csv")

deg1_up = deg1 |> filter(abs(logFC) >1 & P.Value <0.05 & logFC>0)|> pull(symbol)
deg1_down = deg1 |> filter(abs(logFC) >1 & P.Value <0.05 & logFC<0) |> pull(symbol)
deg2_up = deg2 |> filter(abs(logFC) >1 & P.Value <0.05 & logFC>0) |> pull(symbol)
deg2_down = deg2 |> filter(abs(logFC) >1 & P.Value <0.05 & logFC<0) |> pull(symbol)

upList = list(deg1_up, deg2_up)
names(upList) = c("GSE10804_up", "GSE25724_up")
downList = list(deg1_down, deg2_down)
names(downList) = c("GSE10804_down", "GSE25724_down")

ggVennDiagram(upList) + scale_fill_gradient(low="grey90",high = "red")
ggVennDiagram(downList) + scale_fill_gradient(low="grey90",high = "red")


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



ggVennDiagram(sets,
              category.names = c("set 1","set 2","set 3","set 4"), # 集合名字
              set_color = c("blue","red","yellow","black"), # 集合名字颜色
              set_size = 6,# 集合名字大小
              label = "both", # "count","percent","both","none"
              label_geom = "label", # label/text,
              label_alpha = 0, # label背景色
              label_color = "firebrick",
              label_percent_digit = 2, # 保留几位小数
              edge_lty = "dashed", # 边框线型，solid
              edge_size = 1.2 # 边框粗细
              )+
  scale_fill_gradient(low = "grey90",high = "grey60")+ # 填充色
  scale_color_manual(values = c("grey10","grey10","grey10","grey10")) # 边框色,貌似不能变为无边框
