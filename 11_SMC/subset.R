library(tidyverse)
library(Seurat)
library(qs2)
library(future)
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存
rm(list = ls())
setwd("./11_SMC/")

seurat_annotation = qs_read("../09_annotation/seurat_annotation.qs2")

# subset
Idents(seurat_annotation) <- seurat_annotation@meta.data$main_cell_type
table(seurat_annotation$main_cell_type)
SMC <- subset(seurat_annotation, idents = "SMC")

# rescale
# normalize是取对数等操作，基于单元格，不用重新做，但是 scale 操作是基于基因的，所以要重新做，但是 SCT 反推出来的东西靠谱吗？
SMC <- ScaleData(SMC)

# PCA
SMC <- RunPCA(SMC)
ElbowPlot(SMC, ndims = 50)

# save
qs_save(SMC, "SMC.qs2")
