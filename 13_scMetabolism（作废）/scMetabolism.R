library(tidyverse)
library(Seurat)
library(qs2)
library(future)
library(viridis)
library(scMetabolism)
library(rsvd)
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存
setwd("./13_scMetabolism/")
SMC = qs_read("../11_SMC/SMC_cluster.qs2")
Idents(SMC) <- "RNA_snn_res.0.1"

source("./sc.metabolism.SeuratV5.R")
source("./DimPlot.metabolismV5.R")

res <-sc.metabolism.SeuratV5(obj = SMC,
                             method = "VISION", # VISION、AUCell、ssgsea和gsva
                             imputation =F, ncores = 5, 
                             metabolism.type = "KEGG") # KEGG和REACTOME


# check一下有哪些pathway
pathways <- res@assays$METABOLISM$score
head(rownames(pathways))
# [1] "Glycolysis / Gluconeogenesis"            
# [2] "Citrate cycle (TCA cycle)"               
# [3] "Pentose phosphate pathway"               
# [4] "Pentose and glucuronate interconversions"
# [5] "Fructose and mannose metabolism"         
# [6] "Galactose metabolism"

library(wesanderson)
# Dimplot
DimPlot.metabolismV5(obj = res, 
                   pathway = "Glycolysis / Gluconeogenesis", 
                   dimention.reduction.type = "tsne",  
                   dimention.reduction.run = F, size = 1)
phenotype = "RNA_snn_res.0.1"
# DotPlot
input.pathway<-c("Glycolysis / Gluconeogenesis",
                 "Oxidative phosphorylation",
                 "Citrate cycle (TCA cycle)")

input.pathway = rownames(pathways)

DotPlot.metabolism(obj = res,
                  pathway = input.pathway,
                  phenotype = phenotype,
                  norm = "y")
ggsave("./DotPlot_scMetabolism.png", dpi=300, height=10, width=20)

# BoxPlot
# 由于开发者默认吧obj定义为countexp.Seurat，所以还需要重新命名一下
countexp.Seurat <- res
BoxPlot.metabolism(obj = countexp.Seurat,
                  pathway = input.pathway, 
                  phenotype = phenotype, #这个参数需按需修改
                  ncol = 1)
