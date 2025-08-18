library(dplyr)
rm(list = ls())
setwd("./06_PreCyto/")

### net.network 网络关系文件
string_inter = read.table("../05_PPI/string_interactions_short.tsv", header = F)

network = string_inter |> dplyr::select(c(V1, V2))
names(network)  = c("Node1", "Node2")
network$PPI = "ppi"

write.table(network, file="net.network.txt", sep="\t", quote=F, row.names = F)



### 节点属性文件

# 获得上调基因和下调基因的列表
deg1 = read.csv("../01_bulk_anno_diff/GSE10804_deg.csv")
deg2 = read.csv("../01_bulk_anno_diff/GSE25724_deg.csv")

deg1_up = deg1 |> filter(abs(logFC) >0.5 & P.Value <0.05 & logFC>0)|> pull(symbol)
deg1_down = deg1 |> filter(abs(logFC) >0.5 & P.Value <0.05 & logFC<0) |> pull(symbol)
deg2_up = deg2 |> filter(abs(logFC) >0.5 & P.Value <0.05 & logFC>0) |> pull(symbol)
deg2_down = deg2 |> filter(abs(logFC) >0.5 & P.Value <0.05 & logFC<0) |> pull(symbol)

# netnode
genes = c(string_inter$V1, string_inter$V2) |> unique()

upList = list(deg1_up, deg2_up, genes)
downList = list(deg1_down, deg2_down, genes)



# 输出交集基因文件
upInterGenes=Reduce(intersect, upList)        #交集上调基因,Reduce(intersect, list(v1, v2, v3))
upTab=cbind(upInterGenes, "up")
downInterGenes=Reduce(intersect, downList)    #交集下调基因
downTab=cbind(downInterGenes, "down")
interTab=rbind(upTab, downTab)
colnames(interTab)=c("Node", "Type") # 合并
#输出交集基因的属性文件
write.table(interTab, "not.node.txt", sep="\t", quote=F, col.names=T, row.names=F)
#输出基因的列表文件
write.table(upTab, "net.upList.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(downTab, "net.downList.txt", sep="\t", quote=F, col.names=F, row.names=F)
