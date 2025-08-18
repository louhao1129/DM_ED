table(scRNA1@active.ident)

deg=FindMarkers(scRNA1,ident.1 = "Adipocytes",ident.2 = "Granulocytes",
                min.pct = 0.01,logfc.threshold = 0.01)

genelist=deg$avg_log2FC
names(genelist)=toupper(rownames(deg))
# Changing to Upper case.
#result <- toupper("Changing To Upper")
#print(result)
# Changing to lower case.
#result <- tolower("Changing To Lower")
genelist=sort(genelist,decreasing = T)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(enrichplot)
geneset=read.gmt("c2.cp.kegg.v7.4.symbols.gmt")
length(unique(geneset$term))
egmt=GSEA(genelist,TERM2GENE = geneset,
          minGSSize = 1,pvalueCutoff = 0.5)

gseaplot2(egmt,geneSetID = 1,pvalue_table = T)

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



#RcisTarget是一个R-package用于识别基因列表中转录因子(TF)的结合基序（binding motifs ）

#BiocManager::install("RcisTarget")
library(RcisTarget)
# 找到调控marker基因的转录因子
# 这里提取了MA、FI这两类细胞的marker（差异）基因，仅作为展示
geneLists <- list('MA'=deg$gene[1:100],'FI'=deg$gene[200:300]) 

geneLists
#motif 注释文件
data("motifAnnotations_mgi") # 人类用hgs、鼠用mgi
motifAnnotations_mgi
motifRankings <- importRankings("./refdata/mm9-tss-centered-10kb-7species.mc9nr.feather")
##一旦加载了基因列表和数据库，cisTarget()就可以使用它们。cisTarget()运行按顺序执行的步骤
##(1)motif 富集分析,
##(2)motif-TF注释,和
##(3)选择的重要基因。
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot=motifAnnotations_mgi)


motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]

library(DT)
datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))


# 也可以将这些步骤作为单独的命令分别运行。对用户感兴趣的一个输出进行分析，
# 或者优化工作流以在多个基因列表上运行它。我们将分步骤来演示。
#motif 富集分析
##估计每个基序在基因集上的过表达（over-representation）的第一步是计算每对基序-基序集的曲线下面积(AUC)。
#这是根据基因集对基序排序的恢复曲线计算的(根据基序在其邻近度上的得分，基因依次递减，
##                     如motifRanking数据库中提供的那样)。

motifs_AUC <- calcAUC(geneLists, motifRankings, nCores=1)
motifs_AUC
##AUC是由geneset提供的图形矩阵。原则上，AUC主要是作为下一步的输入。
#然而，也有可能探索分数的分布，例如:
par(mfrow = c(1,2))

for(i in c("MA",'FI')){
  auc <- getAUC(motifs_AUC)[i,]
  hist(auc, main=i, xlab="AUC histogram",
       breaks=100, col="#ff000050", border="darkred")
  nes3 <- (3*sd(auc)) + mean(auc)
  abline(v=nes3, col="red")
}


###motif-TF注释

# 显著性motif的选择是基于归一化富集评分( Normalized Enrichment Score,NES)进行的。
# 每个motif的NES是根据基因集中所有基序的AUC分布[(x-mean)/sd]计算。
# 那些通过给定阈值(默认为3.0)的motifs 被认为是重要的。
# motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
#                                            motifAnnot=motifAnnotations_mgi )
 
# ##找出每个基序富集最佳的基因
# 从RcisTarget搜索motif 的富集结果，找到一个motif 是“enriched”并不意味着所有的gene-list的motif 高分。
# 这样，工作流程的第三步是确定哪些基因(基因集)对每个重要的motifs是高度排序的。
# 有两种方法来识别这些基因:(1)相当于那些用于iRegulon和i-cisTarget(method=“iCisTarget”方法,
#  如果运行时间不是问题，建议使用该方法),
# 和(2)更快的一种是实现基于一个近似使用平均分布在每一个等级(method=“aprox”方法,扫描多个基因集)。
 
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                   rankings=motifRankings, 
                                                   geneSets=geneLists)



write.csv(motifEnrichmentTable_wGenes,file = "auc.csv")


geneSetName <- "MA"
selectedMotifs <-sample(motifEnrichmentTable$motif, 3)
par(mfrow=c(2,2))
getSignificantGenes(geneLists[[geneSetName]], 
                    motifRankings,
                    signifRankingNames=selectedMotifs,
                    plotCurve=TRUE, maxRank=5000, genesFormat="none",
                    method="aprox")

# 图中红线为各motif的恢复曲线平均值，绿线为平均值+标准差，蓝线为当前motif的恢复曲线。
# 当前motif与绿色曲线的最大距离点(mean+sd)，为选择的最大富集等级。
# 所有等级较低的基因都被认为是富集的。

#RcisTarget的最终输出是一个包含以下领域中基序富集及其注释信息的表:
motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)

resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]
##可视化
library(DT)
datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))
 
