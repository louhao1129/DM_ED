rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggpubr)


pvalueFilter=0.05     
qvalueFilter=1         

colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd("./03_GO/")       
rt=read.table("../02_venn/interGenes.List.txt", header=F, sep="\t", check.names=F)    

# 获取交集基因名称，并转换为id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       # 去除id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

# 代理
Sys.setenv(http_proxy = "http://127.0.0.1:7890") 
Sys.setenv(https_proxy = "http://127.0.0.1:7890")
#GO
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
saveRDS(kk, "kk.rds")
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)
# 柱状图
pdf(file="barplot.pdf", width=12, height=7)
bar=barplot(kk, drop=TRUE, showCategory=10, label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		
#气泡图
pdf(file="bubble.pdf", width=12, height=7)
bub=dotplot(kk, showCategory=10, orderBy="GeneRatio", label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()
