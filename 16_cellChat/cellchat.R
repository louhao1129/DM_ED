library(CellChat) # https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html#visualize-each-signaling-pathway-using-hierarchy-plot-circle-plot-or-chord-diagram
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)
library(qs2)
library(patchwork)
rm(list = ls())

setwd("./16_cellChat")
SMC = qs_read("../11_SMC/SMC_cluster.qs2")

SMC@meta.data <- SMC@meta.data |> 
  mutate(SMC_cluster = paste0("cluster_", RNA_snn_res.0.1)) |> 
  mutate(SMC_cluster = as_factor(SMC_cluster))

SMC@meta.data$samples = SMC@meta.data$orig.ident
Idents(SMC) = "SMC_cluster"
levels(SMC)
# SMC <- RunUMAP(SMC, dims = 1:20)
# DimPlot(SMC, reduction = "umap", group.by = "RNA_snn_res.0.1")

# 创建cellchat对象
data.input <- LayerData(SMC, assay = "RNA", layer = "data") # normalized data matrix
labels <- Idents(SMC) # 分组信息决定后续cellchat分析的是谁和谁的互作
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = SMC, group.by = "ident", assay = "RNA")

# 设置配体-受体相互作用数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB) # 不出图的时候运行 dev.new()
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)

# 在CellChat中，我们还可以先择特定的信息描述细胞间的相互作用，
# 可以理解为从特定的侧面来刻画细胞间相互作用，比用一个大的配体库又精细了许多。
# 查看可以选择的侧面
unique(CellChatDB$interaction$annotation)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)

cellchat@DB <- CellChatDB.use # set the used database in the object

#对表达数据进行预处理
##将信号基因的表达数据进行子集化，以节省计算成本
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers =5)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
# (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- smoothData(cellchat, adj = PPI.human)

# 细胞间通讯网络推理
## 1、计算通信概率推断细胞互作的通信网络
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
##如果特定细胞群中只有少数细胞，则过滤掉细胞间的通信
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 提取推断出的细胞互作的通信网络数据框
# 函数 Communication，以轻松访问推断出的感兴趣的细胞-细胞通信。例如
##返回一个数据框，包含所有推断的配体/受体级别的细胞-细胞通信。设置slot.name = "netP"以访问信令路径级别的推断通信
df.net <- subsetCommunication(cellchat)
## 给出从 1 和 2 到细胞组 4 和 5 发送的推断细胞间通信
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
## 给出由信号转导 WNT 介导的推断细胞间通讯 和 TGFb
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))

##2、在信号通路水平上推断细胞间的通讯
cellchat <- computeCommunProbPathway(cellchat)
##汇总通信概率来计算细胞间的聚合通信网络。
cellchat <- aggregateNet(cellchat)
qs_save(cellchat, "SMC_cellchat.qs2")

cellchat = qs_read("./SMC_cellchat.qs2")
##3、计算聚合细胞互作通信网络
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#左图：外周各种颜色圆圈的大小表示细胞的数量，圈越大，细胞数越多。发出箭头的细胞表达配体，
#箭头指向的细胞表达受体。配体-受体对越多，线越粗。
#右图：互作的概率或者强度值（强度就是概率值相加）


##每个细胞如何跟别的细胞互作（互作的强度或概率图）
mat <- cellchat@net$weight
par(mfrow = c(1,3), xpd=TRUE) # 根据细胞群数量设置画板
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
##每个细胞如何跟别的细胞互作（number+of+interaction图）
mat <- cellchat@net$count
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


# 可视化：整体、单个通路、配体对

## 使用层次图、圆图或弦图可视化每个信号通路
##查看通路

levels(cellchat@idents)            #查看细胞顺序
vertex.receiver = c(1, 2)          #指定靶细胞的索引
cellchat@netP$pathways             #查看富集到的信号通路
pathways.show <- c("MK")            #指定需要展示的通路

##层次图
vertex.receiver = seq(1,2) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
# 在层次图中，实体圆和空心圆分别表示源和目标。圆的大小与每个细胞组的细胞数成比例。线越粗，互作信号越强。
# 左图中间的target是我们选定的靶细胞。右图是选中的靶细胞之外的另外一组放在中间看互作。

##圈图
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

##和弦图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", vertex.size = groupSize)

##热图
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
##纵轴是发出信号的细胞，横轴是接收信号的细胞，热图颜色深浅代表信号强度。
##上侧和右侧的柱子是纵轴和横轴强度的累积



# 配体-受体层级的可视化（计算各个ligand-receptor+pair对信号通路的贡献）
netAnalysis_contribution(cellchat, signaling = pathways.show)
##也可以看到单个配体-受体对介导的细胞-细胞通信。
#还可以可视化由单配体受体对:我们提供了extractEnrichedLR来提取给定信号通路的所有重要相互作用(L-R对)和相关信号基因。
pairLR.MK <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
#提取对这个通路贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
LR.show <- pairLR.MK[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,2) # a numeric vector
##层次图
netVisual_individual(cellchat, signaling = "MK",  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout="hierarchy")
##圈图
netVisual_individual(cellchat, signaling = "MK", pairLR.use = LR.show, layout = "circle")
##和弦图
netVisual_individual(cellchat, signaling ="MK", pairLR.use = LR.show, layout = "chord")


##############批量保存
pathway.show.all=cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver=c(1,2,3,4)
dir.create("D:/shangke/lession18/cellchat.res")
setwd("D:/shangke/lession18/cellchat.res/")
for (i in 1:length(pathway.show.all)) {

  netVisual(cellchat,signaling = pathway.show.all[i],out.format = c("pdf"),
  vertex.receiver=vertex.receiver,layout="circle")
plot=netAnalysis_contribution(cellchat,signaling = pathway.show.all[i])
  ggsave(filename = paste0(pathway.show.all[i],".contribution.pdf"),
         plot=plot,width=6,height=4,dpi=300,units="in")
  
}

#####################################################################################
# 可视化由多个配体受体或信号通路介导的细胞间通讯
##气泡图
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:3), remove.isolate = FALSE)
##sources.use = 2 是值第二个细胞亚群
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use = c(1:3), remove.isolate = FALSE)
##指定信号通路
cellchat@netP$pathways 
netVisual_bubble(cellchat, sources.use =c(1,2), targets.use =c(1:3),signaling =  c("MK","CypA"), remove.isolate = FALSE)


pairLR  <- extractEnrichedLR(cellchat, signaling =c("MK","CypA"), geneLR.return = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,2), targets.use =c(1:3),pairLR.use =pairLR , remove.isolate = FALSE)

#和弦图               
netVisual_chord_gene(cellchat, sources.use = 2, targets.use = c(1:3), lab.cex = 0.5,legend.pos.y = 30)
# 更多见官方文档
##用小提琴图绘制信号基因的表达分布 参与某条信号通路（如TGFb）的所有基因在细胞群中的表达情况展示
plotGeneExpression(cellchat, signaling = "MK")
#默认情况下，plotGeneExpression只显示与推断的重要通信相关的信号基因的表达。
#用户可以通过显示一个信号通路相关的所有信号基因的表达。
plotGeneExpression(cellchat, signaling = "MK", enriched.only = FALSE) # 注意enriched.only参数
#也可以用气泡图展示
plotGeneExpression(cellchat, signaling = "MK",type = "dot")

###################################################################################################
# 细胞通信网络系统分析
##可视化配体和受体
## 1、计算网络中心性得分
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
##2、热图  使用热图可视化计算的中心性评分，允许随时识别细胞群的主要信号作用。
netAnalysis_signalingRole_network(cellchat, signaling = "MK", width = 8, height = 2.5, font.size = 10)

##在2D空间中可视化主要的发送者(源)和接收者(目标)。
##我们还提供了另一种直观的方式，使用散点图来可视化2D空间中的主要发送者(源)和接收者(目标)。

##从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
###从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MK", "PARs"))
gg1 + gg2

##识别对某些细胞群的传出或传入信号贡献最大的信号，从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析。

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))
########################################################################################################
# 识别全局通讯模式，探索多种细胞类型和信号通路如何协调在一起
# 细胞通讯模式和信号网络
# 除了探索单个通路的详细通信外，一个重要的问题是多个细胞群和信号通路如何协调发挥作用。CellChat 采用模式识别方法来识别全局通信模式。
library(NMF)
library(ggalluvial)
#非负矩阵分解（NMF）识别细胞的通讯模式
##信号输出细胞的模式识别
##计算分解成几个因子(pattern)比较合适（这一步运行比较慢+。在使用NMF对细胞进行亚群细分时，如果不测试的话，最好选择比细胞类型多一点的值）
selectK(cellchat, pattern = "outgoing")


#挑选曲线中第一个出现下降的点（从3就开始下降了）
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
##river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#气泡图
netAnalysis_dot(cellchat, pattern = "outgoing")
#信号输入细胞的模式识别
selectK(cellchat, pattern = "incoming")

#################################################################################################
#  信号网络聚类
# 1、根据功能相似性来识别信号分组

##reticulate::py_install(packages = 'umap-learn')
 
## 2、基于结构相似性识别信号分组
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

#############################################################################################
#不同分组之间的配对分析
# 疾病组合对照组的配对分析，创建两个cellchat对象，分别分析，然后merge，进行可视化
sc.sp=SplitObject(scRNA_harmony,split.by = "orig.ident")
sc.11=scRNA_harmony[,sample(colnames(sc.sp[["sample_11"]]),1000)]
sc.3=scRNA_harmony[,sample(colnames(sc.sp[["sample_3"]]),1000)]


 
cellchat.sc11 <- createCellChat(object =sc.11@assays$RNA@data, meta =sc.11@meta.data,  group.by ="celltype")
cellchat.sc3 <- createCellChat(object =sc.3@assays$RNA@data, meta =sc.3@meta.data,  group.by ="celltype")

dir.create("compare")
setwd("compare/")

cellchat=cellchat.sc11 
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc11 = cellchat
#################################
cellchat=cellchat.sc3
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc3 = cellchat
##############################################
cc.list=list(SC11=cc.sc11,SC3=cc.sc3)
cellchat=mergeCellChat(cc.list,cell.prefix = T,add.names = names(cc.list))
##可视化
##所有细胞群总体观：通讯数量与强度对比
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "count")
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "weight")
##第一个图展示通讯数量之间的差异，第二个图展示通讯强度之间的差异。 

##数量与强度差异网络图
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight")
##红色是case相对于control上调的，蓝色是下调的

#数量与强度差异热图
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat,measure = "weight")
#case和control对比，红色是上调，蓝色是下调

#保守和特异性信号通路的识别与可视化
rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
rankNet(cellchat,mode = "comparison",stacked =F,do.stat = T)
##左图最下面多个信号通路是case组独有的

##细胞互作数量对比网络图
weight.max=getMaxWeight(cc.list,attribute = c("idents","count"))
netVisual_circle(cc.list[[1]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )

netVisual_circle(cc.list[[2]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )


table(scRNA_harmony@active.ident)
s.cell=c( "Macrophage", "Tissue_stem_cells","Monocyte")
count1=cc.list[[1]]@net$count[s.cell,s.cell]
count2=cc.list[[2]]@net$count[s.cell,s.cell]

netVisual_circle(count1,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )

netVisual_circle(count2,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )

