# pyscenic参考资料
# 生信技能树：https://cloud.tencent.com/developer/article/2228252
# Biomamba：https://zhuanlan.zhihu.com/p/31567655578?share_code=ACwaAxFVYgGJ&utm_psn=1933240631137920295
# cisTarget 资源网站：https://resources.aertslab.org/cistarget/

# SCENIC/pySCENIC 评估转录因子活性，其实也是通过这个转录因子调控的基因的表达来评估的，
# 可以理解为评估一条通路（这里是一个转录因子）的活性，类似基因集富集分析(GSEA)，不过使用的算法是AUCell

##可视化
rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
setwd("./14_pyscenic")

#### 1.提取 out_SCENIC.loom 信息
loom <- open_loom('results/out_SCENIC.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons") #获取调节因子矩阵
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat) #将调节因子矩阵转换为基因列表
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC') #获取调节因子的AUC值
regulonAucThresholds <- get_regulon_thresholds(loom) #获取调节因子的阈值
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)

#### 2.加载SeuratData
SMC = qs_read("../11_SMC/SMC_cluster.qs2")
Idents(SMC) = "RNA_snn_res.0.1"
levels(SMC)
DimPlot(SMC, reduction="tsne", label = T, label.size = 5, split.by = "group")

seurat.data = SMC
sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]
dim(sub_regulonAUC)
seurat.data
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data$RNA_snn_res.0.1))
cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$RNA_snn_res.0.1)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

#保存一下
save(sub_regulonAUC,cellTypes,cellClusters,seurat.data,
     file = 'for_rss_and_visual.Rdata')

load("./for_rss_and_visual.Rdata")


# #每个细胞的aucell值添加到seurat对象的meta.data中
seurat.data@meta.data = cbind(seurat.data@meta.data ,t(assay(sub_regulonAUC)))

# 指定关注的转录因子
regulonsToPlot = c('NR2F1(+)','PITX1(+)','FOXP1(+)')
regulonsToPlot %in% row.names(sub_regulonAUC)

# Vis
p1 = DotPlot(seurat.data, features = unique(regulonsToPlot)) + RotatedAxis()
p2 = RidgePlot(seurat.data, features = regulonsToPlot , ncol = 2) 
p3 = VlnPlot(seurat.data, features = regulonsToPlot,pt.size = 0)
p4 = FeaturePlot(seurat.data,features = regulonsToPlot, reduction = 'tsne')

wrap_plots(p1,p2,p3,p4)
ggsave("./regulonsToPlot.png", dpi=300, width=15, height=10)
#可选
#scenic_res = assay(sub_regulonAUC) %>% as.matrix()
#seurat.data[["scenic"]] <- SeuratObject::CreateAssayObject(counts = scenic_res)
#seurat.data <- SeuratObject::SetAssayData(seurat.data, slot = "scale.data",
#                                  new.data = scenic_res, assay = "scenic")

### 4.1. TF活性均值
# 看看不同单细胞亚群的转录因子活性平均值
# Split the cells by cluster:
selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

# 去除extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression. 
# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
# 参考：https://www.jianshu.com/p/115d07af3029
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)
#[1] 209   9
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

png("regulon_heatmap.png", width = 1500, height = 5000, res = 300)
Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 5),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

dev.off()
# 可以看到，确实每个单细胞亚群都是有自己的特异性的激活的转录因子。

### 4.2. rss查看特异TF
# 不过，SCENIC包自己提供了一个 calcRSS函数，帮助我们来挑选各个单细胞亚群特异性的转录因子，全称是：Calculates the regulon specificity score
# regulon特异性分数(Regulon Specificity Score, RSS)
# 参考文章：The RSS was first used by Suo et al. in: Revealing the Critical Regulators of Cell Identity in the Mouse Cell Atlas. Cell Reports (2018). doi: 10.1016/j.celrep.2018.10.045 运行超级简单。
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)
rssPlot$plot
plotly::ggplotly(rssPlot$plot)

# 分细胞类型展示regulon特异性分数
plotRSS_oneSet(rss,'0')
plotRSS_oneSet(rss,'1')
plotRSS_oneSet(rss,'2')
# 4.3 其他查看TF的方式
# 从 regulon 活性打分结果中，找出每个 cluster 中相对于其他 cluster 活性最显著的前 5 个 TF regulon
# 并可视化这些高特异性的 TF 活性热图
rss=regulonActivity_byGroup_Scaled
head(rss)
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),           # TF regulon 名
                 cluster = colnames(rss)[i],      # 当前 cluster 名
                 sd.1 = rss[,i],                  # 当前 TF 在当前 cluster 的得分
                 sd.2 = apply(rss[,-i], 1, median)  # 当前 TF 在其他 cluster 的中位数得分
               )
             }))
df$fc = df$sd.1 - df$sd.2 # 计算 TF 的“差异得分”（fold change 概念）
# 这里的 fc 是每个 regulon 在该 cluster 中的特异性打分，高说明特异性强。
#挑选每个 cluster 的 top 5 高特异性 TF regulon
top5 <- df %>% group_by(cluster) %>% top_n(5, fc) 

rowcn = data.frame(path = top5$cluster)  # 行注释
n = rss[top5$path,]                      # 提取 top TF 的 regulon 活性

pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T)
