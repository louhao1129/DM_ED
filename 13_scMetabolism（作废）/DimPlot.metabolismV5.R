# 需要修改DimPlot.metabolism中的UMAP大小写
DimPlot.metabolismV5 <- function (obj, pathway, dimention.reduction.type = "umap", dimention.reduction.run = T, 
  size = 1) 
{
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  if (dimention.reduction.type == "umap") {
    if (dimention.reduction.run == T) 
      obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:40)
    umap.loc <- obj@reductions$umap@cell.embeddings
    row.names(umap.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(umap.loc, t(signature_exp[input.pathway, 
      ]))
    library(wesanderson)
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot <- ggplot(data = signature_ggplot, aes(x = umap_1, 
      y = umap_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal) + 
      labs(color = input.pathway) + xlab("UMAP 1") + ylab("UMAP 2") + 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), panel.background = element_blank(), 
      axis.line = element_line(colour = "black"))
  }
  if (dimention.reduction.type == "tsne") {
    if (dimention.reduction.run == T) 
      obj <- Seurat::RunTSNE(obj, reduction = "pca", dims = 1:40)
    tsne.loc <- obj@reductions$tsne@cell.embeddings
    row.names(tsne.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(tsne.loc, t(signature_exp[input.pathway, 
      ]))
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot <- ggplot(data = signature_ggplot, aes(x = tSNE_1, 
      y = tSNE_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal) + 
      labs(color = input.pathway) + xlab("tSNE 1") + ylab("tSNE 2") + 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), panel.background = element_blank(), 
      axis.line = element_line(colour = "black"))
  }
  plot
}
