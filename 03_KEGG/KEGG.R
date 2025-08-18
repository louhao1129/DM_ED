rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05      
qvalueFilter=1       # 矫正后p值的过滤条件，设置为1说明不按它进行过滤，可设置为0.05 

# 自动判断结果显示p value还是q value
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue" 
}
	
setwd("./03_KEGG/")     
rt=read.table("../02_venn/interGenes.List.txt", header=F, sep="\t", check.names=F)     

# 把基因名转换为id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#KEGG
Sys.setenv(http_proxy = "http://127.0.0.1:7890") 
Sys.setenv(https_proxy = "http://127.0.0.1:7890")

kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
saveRDS(kk, "kk.rds")
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)


showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}


pdf(file="barplot.pdf", width=12, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()


pdf(file="bubble.pdf", width=12, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()

