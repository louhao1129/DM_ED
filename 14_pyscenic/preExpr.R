library(Seurat)
library(tidyverse)
library(patchwork)
library(qs2)

setwd("./14_pyscenic/")
SMC = qs_read("../11_SMC/SMC_cluster.qs2")
Idents(SMC) = "RNA_snn_res.0.1"
levels(SMC)

# 提取表达矩阵
counts <- LayerData(SMC, assay = "RNA", layer = "counts") |> as.matrix()
head(colnames(counts))
head(rownames(counts))
write.csv(counts,file = "./result/counts.csv")
