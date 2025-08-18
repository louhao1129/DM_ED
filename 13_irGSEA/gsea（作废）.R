library(tidyverse)
library(Seurat)
library(qs2)
library(future)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(enrichplot)
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存
rm(list = ls())
setwd("./13_gsea/")
SMC = qs_read("../11_SMC/SMC_cluster.qs2")
# Find Marker
Idents(SMC) = "RNA_snn_res.0.1"
levels(SMC)
# find marker for some cluster
# cluster.markers <- FindMarkers(SMC, ident.1 = 0)
# cluster.markers |> 
#   filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5) & pct.1 >0.7 & pct.2<0.3)

# find markers for every cluster compared to all remaining cells, report only the positive
# SMC.markers <- FindAllMarkers(SMC,logfc.threshold = -Inf, min.pct = 0.01)
# GSEA则不局限于差异基因，不需要也不应该根据 fold change 筛选差异基因。
# GSEA 核心思想：对所有基因按 fold change（或其他统计量）排序，然后判断某个通路的基因是否倾向集中在排序的上端或下端。
marker_cluster_0 = SMC.markers |> 
  filter(cluster==0)
marker_cluster_1 = SMC.markers |> 
  filter(cluster==1)
marker_cluster_2 = SMC.markers |> 
  filter(cluster==2)
marker_cluster_1_2 = SMC.markers |> 
  filter(cluster==2 | cluster==1)

# GSEA只需要基因名和foldchage
marker_cluster = marker_cluster_1_2
genelist=marker_cluster$avg_log2FC
names(genelist) = rownames(marker_cluster)
genelist=sort(genelist,decreasing = T)

geneset=read.gmt("./KEGG_metabolism_nc.gmt") # 其实就是一个dataframe
length(unique(geneset$term)) # 查看通路数量
table(geneset$term) # 查看每个通路的基因数量
egmt=GSEA(
          genelist,
          TERM2GENE = geneset,
          minGSSize = 1, # 每个基因集（gene set）中最少要包含的基因数目，低于这个值的通路将被忽略。
          pvalueCutoff = 0.05)

gseaplot2(egmt,geneSetID = 1,pvalue_table = T)
gseaplot(egmt,geneSetID = 1,pvalue_table = T)
kegg.res=egmt@result
down.kegg.res<-kegg.res[(kegg.res$pvalue<0.05 & kegg.res$enrichmentScore < -0.3),]
down.kegg.res$group=1

up.kegg.res<-kegg.res[(kegg.res$pvalue<0.05 & kegg.res$enrichmentScore > 0.3),]
up.kegg.res$group=1

lapply(1:nrow(down.kegg.res), function(i){
  
  gseaplot2(egmt,down.kegg.res$ID[i],title = down.kegg.res$Description[i],pvalue_table = T)
  ggsave(paste0(gsub("/","-",down.kegg.res$Description[i]),".pdf"),width = 11,height =5)
}
  )


lapply(1:nrow(up.kegg.res), function(i){
  
  gseaplot2(egmt,up.kegg.res$ID[i],title = up.kegg.res$Description[i],pvalue_table = T)
  ggsave(paste0(gsub("/","-",up.kegg.res$Description[i]),".up.pdf"),width = 11,height =5)
}
)
#########################################################################################
##install.packages("msigdbr")
 
library(msigdbr)
##BiocManager::install("fgsea")
library(fgsea)
msigdbr_species()

# fgsea能够快速对预选基因集进行GSEA富集分析，预选基因集可以是自己设定，一般使用MSigdb数据库（同样由提出GSEA方法团队提供）。
# fgsea（）函数需要一个基因集列表以及对应值，主要是基因名和AUC(ROC曲线下方的面积大小，简单说就是代表准确性，准确性越高，AUC值越大)，其中基因集中的基因名要与数据集（pathway）中的基因名相对应。
# fgsea包中的plotEnrichment函数用于GSEA图的绘制。
 
m_df<- msigdbr(species = "Mus musculus",  category = "C2", subcategory = "KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)



genelist=deg$avg_log2FC
names(genelist)= rownames(deg)
 
genelist=sort(genelist,decreasing = T)
?fgsea
fgseaRes<- fgsea(fgsea_sets, stats = genelist )

ggplot(fgseaResTidy %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 1.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal() ####以1.5进行绘图填色


plotEnrichment(fgsea_sets[["KEGG_ABC_TRANSPORTERS"]],
               genelist) + labs(title="KEGG_ABC_TRANSPORTERS")

