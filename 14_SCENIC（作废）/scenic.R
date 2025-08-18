##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC) # https://scenic.aertslab.org/tutorials/
library(harmony)
library(qs2)
setwd("./16_SCENIC/")
SMC = qs_read("../11_SMC/SMC_cluster.qs2")
Idents(SMC) = "RNA_snn_res.0.1"
levels(SMC)

exprMat <- LayerData(SMC, assay = "RNA", layer = "counts")
exprMat <- as.matrix(exprMat)
# 基因 ID 应该是gene symbol并存储为行名 （以便与 RcisTarget 注释数据库兼容）。
dim(exprMat)
head(rownames(exprMat))
##设置分析环境
mydbDIR <- "./scenic"
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")
#初始化 SCENIC 设置,设置分析环境
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

scenicOptions <- initializeScenic(org="hgnc", # 人类
                                  nCores=8,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "os")
saveRDS(scenicOptions, "int/scenicOptions.rds")


##==转录调控网络推断==##
##基因过滤
#过滤标准是基因表达量之和>细胞数*3%，且在1%的细胞中表达
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)

# 是否把我们感兴趣的基因给过滤了？
# interestingGenes <- c("Sox9", "Sox10", "Dlx5")
# interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)
##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions) # 输出在int文件夹中
##TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
# 根据表达数据推断潜在的转录因子靶标，使用 GENIE3 或 GRNBoost，
# GENIE3 非常耗时且计算量大（在 3-5k 单元的数据集上需要几个小时或几天的时间）
# GRNboost可在很短的时间内提供与 GENIE3 类似的结果，这儿使用的R，选择GENIC3

##nParts参数，是把表达矩阵分成n份分开计算，防止内存不够
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
head("int/1.4_GENIE3_linkList.Rds") # 这是GENIE3最终的结果，是把1.3_开头的文件（也就是中间结果）合并到了一起
head(readRDS("int/1.4_GENIE3_linkList.Rds") )
# TF是转录因子的名称，Target是潜在靶基因的名称，weight是TF和Target之间相关性的权重

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) #1. 获取共表达模块,保存为1.6_tfModules_asDF.Rds
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)  #2. 获取regulons（转录调控单元）
# 2.1对每个共表达模块进行motif富集分析，保留显著富集的motif；此项分析依赖gene-motif评分（排行）数据库，其行为基因、列为motif、值为排名，就是我们下载的cisTarget数据库。
# 2.2使用数据库对motif进行TF注释，注释结果分高、低可信度。数据库直接注释和同源基因推断的TF是高可信结果，使用motif序列相似性注释的TF是低可信结果。
# 2.3用保留的motif对共表达模块内的基因进行打分（同样依据cisTarget数据库），识别显著高分的基因（理解为motif离这些基因的TSS很近）；
# 2.4删除共表达模块内与motif评分不高的基因，剩下的基因集作者称为调控单元。
# TSS搜索空间:10kbaroundtheTSSor50obpupstreamtheTSS
# 保存的结果在output文件夹
?runSCENIC_2_createRegulons

##==regulon活性评分与可视化==##
##regulons计算AUC值并进行下游分析
load("scRNA_harmony.Rdata")
mydbDIR <- "D:/ref.data/"
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")
#初始化 SCENIC 设置,设置分析环境
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "os")
library(foreach)
exprMat_all <- as.matrix(scRNA_harmony@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
rm(scRNA_harmony)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)
# runSCENIC_3_scoreCellsO是一个多种功能打包的函数，它包含有三个子功能，分别为：
# 3.1计算regulon在每个细胞中AUc值，最后得到一个以regulon为行细胞为列的矩阵。
# 结果目录：int/3.4_regulonAUC.Rds
# 3.2计算每个regulon的AUC阈值，细胞中regulonAUC值>阈值，代表此regulon在细胞中处于激活状态，否则代表沉默状态。
# 结果目录：5_AUCellThresholds.Rds
# 3.3使用regulonAUC矩阵对细胞进行降维聚类
# 用heatmap图展示regulonAUC矩阵，用t-SNE图分别展示每个regulon的活性分布情况。
# 结果目录：output/Step3_开头的系列文件

#使用shiny互动调整阅值
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_all)
# savedSelections <-shiny::runApp(aucellApp)
# #保存调整后的阈值
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1]
# "int/newThresholds.Rds"
# saveRDS(newThresholds,file=getIntName(scenicOptions,"aucell_thresholds"))
# saveRDS(scenicoptions, file="int/scenicOptions.Rds")

runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all) 
# 4、Regulon活性二进制转换与可视化
# 对于细胞类型清晰的数据集而言，构建regulons并对其活性打分足够后续分析。
# 但是，在很多情况下将评分转换为二进制的"开/关”，则既方便解释，又能最大化体现细胞类型的差异。
# 将特定的regulon转换为"0/1"有利于探索和解释关键转录因子。将所有的regulons转换为"0/1"后创建二进制的活性矩阵，则可以用于细胞聚类，对消除技术偏倚特别有用。
# AUCel会自动计算可能的阀值进行"0/1"转换，作者建议在转换之前手工检查并调整这些阈值。
# 将regulonAUC矩阵转换为二进制矩阵后，会重新降维聚类，输出的结果与runSCENIC3scoreCells（）类似

##准备细胞meta信息，一定要按下面的格式保存到int文件夹下，方便SENIC识别
cellInfo <- data.frame(scRNA_harmony@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype_Monaco")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster","celltype")]
saveRDS(cellInfo, file="int/cellInfo.Rds")

mydbDIR <- "D:/ref.data/"
# 着重看样本的处理方式，使用了哪个参考基因组
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")
#初始化 SCENIC 设置,设置分析环境
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "os")


library(foreach)
nPcs <- c(5)

fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))

# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="celltype",   cex=.5)



 
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")


library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_all, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("JUN","MYC")],], plots="Expression")

#par(bg = "black")
par(mfrow=c(1,2))
regulonNames <- c( "JUN","MYC")
cellCol <- plotEmb_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)


regulonNames <- list( green=c("JUN"),
                     blue=c( "MYC"))
cellCol <- plotEmb_rgb(scenicOptions, regulonNames, aucType="Binary")

regulons <- loadInt(scenicOptions, "regulons")
regulons[c( "JUN","MYC")]

regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="YY1" & highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=5)) 

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="YY1"]
viewMotifs(tableSubset)



regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
 
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:20,], name="Regulon activity")


topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$celltype), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
ComplexHeatmap::Heatmap(binaryActPerc_subset[1:20,], name="Regulon activity (%)", col = c("white","pink","red"))


topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
viewTable(topRegulators)



rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "celltype"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)


plotRSS_oneSet(rss, setName = "MSC")


library(Seurat)
#scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:16)
dr_coords <- Embeddings(scRNA_harmony, reduction="tsne")
par(mfrow=c(4,4))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, "MYC"), plots = "AUC")


tfs <- c("HDAC2","RAD21","YY1", "SMARCA4")
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")
########################################################################


##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(scRNA_harmony, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA_harmony, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'scRNAbin.rds')

##利用Seurat可视化AUC
dir.create('scenic_seurat')
#FeaturePlot
colnames(scRNAauc@meta.data)[20:30]
p1 = FeaturePlot(scRNAauc, features="HCFC1_24g", label=T, reduction = 'tsne')
p2 = FeaturePlot(scRNAbin, features="HCFC1_24g", label=T, reduction = 'tsne')
p3 = DimPlot(scRNA_harmony, reduction = 'tsne', group.by = "celltype", label=T)
plotc = p1|p2|p3
plotc



#RidgePlot&VlnPlot
p1 = RidgePlot(scRNAauc, features = "NR3C1_1339g", group.by="celltype") + 
        theme(legend.position='none')
p2 = VlnPlot(scRNAauc, features ="NR3C1_1339g", pt.size = 0, group.by="celltype") + 
        theme(legend.position='none')
plotc = p1 + p2
plotc


library(pheatmap)
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select = 'celltype')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
#挑选部分感兴趣的regulons
my.regulons <- rownames(AUCmatrix)[80:100]
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]
#使用regulon原始AUC值绘制热图
pheatmap(myAUCmatrix, show_colnames=F, annotation_col=celltype )
#使用regulon二进制AUC值绘制热图
pheatmap(myBINmatrix, show_colnames=F, annotation_col=celltype,
         color = colorRampPalette(colors = c("white","black"))(100))