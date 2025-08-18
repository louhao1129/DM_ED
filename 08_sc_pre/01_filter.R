rm(list = ls())
library(data.table)
library(Seurat)
library(qs2)
library(future)
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存
setwd("./08_sc_pre/")

dir="./GSE206528_RAW/"
samples=list.files( dir )
samples

sceList = lapply(samples,function(pro){ 
  print(pro) # 打印当前处理的文件名
  ct=fread(file.path(dir ,pro),data.table = F)
  ct[1:4,1:4] # 预览数据前4行4列（调试用）
  rownames(ct)=ct[,1] # 将第一列设为行名（基因名）
  ct=ct[,-1] # 删除已作为行名的第一列
  sce=CreateSeuratObject(counts =  ct ,
                         project = gsub('_gene_matrix.csv.gz','',gsub('^GSM[0-9]*_','',pro) ) ,
                         min.cells = 3, # 保留在>=3个细胞中表达的基因
                         min.features = 200) # 保留检测到>=200个基因的细胞
  
  return(sce)
})
samples
#合并数据
sce.all=merge(x=sceList[[1]],y=sceList[ -1 ],
              add.cell.ids =  gsub('_gene_matrix.csv.gz','',gsub('^GSM[0-9]*_','',samples)))
dim(sce.all)
# [1] 22868 64993
# 3 名勃起正常男性和 5 名器质性 ED 患者的 64,993 个单独的海绵细胞

# quality control
nFeature_lower <- 800
nFeature_upper <- 7000
nCount_lower <- 0
nCount_upper <- 30000
pMT_upper <- 10  # 线粒体比例上限

# 计算线粒体基因比例
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")

# 应用所有过滤条件
seurat_filtered <- subset(
  sce.all,
  subset = nFeature_RNA > nFeature_lower & 
           nFeature_RNA < nFeature_upper &
           nCount_RNA > nCount_lower & 
           nCount_RNA < nCount_upper &
           percent.mt < pMT_upper
)
# 22868 64992
qs_save(seurat_filtered, "./seurat_filtered.qs2")