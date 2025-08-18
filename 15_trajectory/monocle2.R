library(monocle)
library(Seurat)
library(tidyverse)
library(patchwork)
library(qs2)
rm(list = ls())

setwd("./15_trajectory/")
SMC = qs_read("../11_SMC/SMC_cluster.qs2")
Idents(SMC) = "RNA_snn_res.0.1"
levels(SMC)

# 提取表达矩阵
counts <- LayerData(SMC, assay = "RNA", layer = "counts") |> as.matrix()
data <- as(counts, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = SMC@meta.data) # 使用new创建S4对象
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.1,
                              expressionFamily = negbinomial.size())



monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))

# 选择高变基因，也可以使用别的方法，选择基因进行轨迹构建
HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)


HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')


HSMM <- orderCells(HSMM, reverse = TRUE)
# 如果 Monocle2 自动推断的拟时序方向与你的生物学认知相反，就设置 reverse = TRUE 来翻转拟时序方向，默认为NULL
# 如果出现igraph版本报错，先重启R终端，然后运行下面的代码安装旧版本
# install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz",repos = NULL)

plot_cell_trajectory(HSMM, color_by = "RNA_snn_res.0.1")
plot_cell_trajectory(HSMM, color_by = "RNA_snn_res.0.1")+
  facet_wrap(~group, nrow = 1)

plot_cell_trajectory(HSMM, color_by = "State") # State是根据monocle算法得到的状态

# 注意这只是算法的结果，Pseudotime如果与生物学结果不符合，要以生物学结果为依据
# monocle3可以自定义发育起点
plot_cell_trajectory(HSMM, color_by = "Pseudotime")

# 查看特定基因的表达量
blast_genes <- row.names(subset(fData(HSMM),
                                gene_short_name %in% c("GAPDH", "RORA")))
plot_genes_jitter(HSMM[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)



HSMM_expressed_genes <-  row.names(subset(fData(HSMM),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("YWHAB", "GAPDH", "TNNC1")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "RNA_snn_res.0.1")
plot_genes_in_pseudotime(cds_subset, color_by =  "State")

genes<-c("TNNT2","TNNC1","CDK1")
p1<-plot_genes_jitter(HSMM[genes,],grouping ="State",color_by="State")
p2<-plot_genes_violin(HSMM[genes,],grouping ="State",color_by="State")
p3 <-plot_genes_in_pseudotime(HSMM[genes,],color_by ="State")
plotc<-p1|p2|p3
plotc

# 根据拟时序表达的模式对基因进行聚类
marker_genes <- row.names(subset(fData(HSMM),
                   gene_short_name %in% c("MEF2C", "MEF2D", "MYF5",
                                          "ANPEP", "PDGFRA","MYOG",
                                          "TPM1",  "TPM2",  "MYH2",
                                          "MYH3",  "NCAM1", "TNNT1",
                                          "TNNT2", "TNNC1", "CDK1",
                                          "CDK2",  "CCNB1", "CCNB2",
                                          "CCND1", "CCNA1", "ID1")))

diff_test_res <- differentialGeneTest(HSMM[marker_genes,],
              fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                num_clusters = 3, # cluster的数量
                cores = 1,
                show_rownames = T)

# 单细胞分支分析

BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 5)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,
                                          qval < 1e-4)),],
                                          branch_point = 1,
                                          num_clusters = 4,
                                          cores = 1,
                                          use_gene_short_name = T,
                                          show_rownames = T)