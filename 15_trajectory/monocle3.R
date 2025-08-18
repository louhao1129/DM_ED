##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle3)
library(patchwork)
library(qs2)

setwd("./15_trajectory/")
SMC = qs_read("../11_SMC/SMC_cluster.qs2")
Idents(SMC) = "RNA_snn_res.0.1"
levels(SMC)
SMC <- RunUMAP(SMC, dims = 1:20)
DimPlot(SMC, reduction = "umap", group.by = "RNA_snn_res.0.1")

data <- LayerData(SMC, assay = "RNA", layer = "counts")
cell_metadata <- SMC@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
#cds <- preprocess_cds(cds, num_dim = 20)
#plot_pc_variance_explained(cds)
reducedDims(cds)$PCA = Embeddings(SMC, "pca") # 使用 Seurat 的 PCA 降维结果

## Step 2(optional): Remove batch effects with cell alignment
#cds <- align_cds(cds, alignment_group = "orig.ident")

## Step 3: Reduce the dimensions 
# cds <- reduce_dimension(cds)
reducedDims(cds)$UMAP = Embeddings(SMC, reduction = "umap") # 使用Seurat的UMAP降维结果

plot_cells(cds, color_cells_by="RNA_snn_res.0.1", group_label_size =5)
plot_cells(cds, color_cells_by="orig.ident", group_label_size =5)
plot_cells(cds, genes=c("GAPDH"))
plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE, group_label_size =5)

cds = cluster_cells(cds) # If NULL (Default), the parameter is determined automatically.
table(cds@clusters$UMAP$clusters)
cds@clusters$UMAP$partitions
# cluster：指的是细胞类型的局部分群，体现细胞之间的表达相似性，和 Seurat 中的 cluster 类似。
# partition：是轨迹结构上的大分区，指的是细胞之间是否共享同一条轨迹树，不同 partition 之间的细胞轨迹是不连通的。

cds@clusters$UMAP$clusters <- factor(SMC$RNA_snn_res.0.1) # 使用Seurat的cluster

plot_cells(cds, color_cells_by="partition")
plot_cells(cds, color_cells_by="cluster", group_label_size =5)

# marker_test_res <- top_markers(cds, group_cells_by="partition", 
#                                reference_cells=1000, cores=8)

# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)

# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="celltype",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)


# top_specific_markers = marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(3, pseudo_R2)

# top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))

# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="celltype",
#                     ordering_type="cluster_row_col",
#                     max.size=3)

# 先将partitions的分组由因子型转为字符型
# colData(cds)$assigned_cell_type <- as.character(partitions(cds))

# cds_subset <- choose_cells(cds)
# 我们已经使用Seurat提取了子集，这里不再提取

cds <- learn_graph(cds, use_partition = F) # learn_graph() 只支持 UMAP，不支持 TSNE
# 这里只有SMC一种类型的细胞，理论上应该共享同一条轨迹树，整张轨迹图应该是是连通的，所有 cluster 都在一张轨迹图上演化
plot_cells(cds,
           color_cells_by = "RNA_snn_res.0.1",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size=4,cell_size=1.5)

cds = order_cells(cds) # only "UMAP" is supported

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)

table(scRNA_harmony@meta.data$celltype)
# a helper function to identify the root principal points:
get_earliest_principal_node  <- function(cds, time_bin="T_cells"){
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)


cds_sub <- choose_graph_segments(cds)

#Working with 3D trajectories

cds_3d = reduce_dimension(cds, max_components = 3)
cds_3d = cluster_cells(cds_3d)
cds_3d = learn_graph(cds_3d)
cds_3d = order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj = plot_cells_3d(cds_3d, color_cells_by="celltype")
cds_3d_plot_obj


ciliated_genes = top_specific_markers$gene_id[5:10]
cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]

gene_fits = fit_models(cds_subset, model_formula_str = "~seurat_clusters")
fit_coefs
fit_coefs = coefficient_table(gene_fits)
# 挑出时间相关的组分
emb_time_terms = fit_coefs %>% filter(term == "seurat_clusters1")
emb_time_terms

emb_time_terms = fit_coefs %>% filter(term == "seurat_clusters1")

emb_time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

plot_genes_violin(cds_subset[,], group_cells_by="celltype", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


neurons_cds <- cds[,grepl("Macrophage", colData(cds)$celltype, ignore.case=TRUE)]
plot_cells(neurons_cds, color_cells_by="partition")

pr_graph_test_res <- graph_test(neurons_cds, neighbor_graph="knn", cores=3)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_module_df = find_gene_modules(neurons_cds[pr_deg_ids,], resolution=1e-2)

cell_group_df = tibble::tibble(cell=row.names(colData(neurons_cds)), cell_group=neurons_cds@colData@listData[["seurat_clusters"]])
agg_mat = aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) = stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)


plot_cells(neurons_cds,
           genes=gene_module_df %>% filter(module %in% c(16,38,33,42)),
           group_cells_by="cluster",
           color_cells_by="cluster",
           show_trajectory_graph=FALSE)


plot_cells(cds,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

 
 ##找到影响发育轨迹的基因
trace('calculateLW', edit = T, where = asNamespace("monocle3"))
## change Matrix::rBind to rbind on line 93
# 使用neighbor_graph="principal_graph"来检验轨迹相邻的细胞的表达是否相关
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)



pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
 
AFD_genes <-pr_deg_ids[20:22]
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                        ]
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="celltype",
                         min_expr=0.5)


cds_subset <- choose_cells(cds)

subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)

agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)




