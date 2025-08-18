rm(list=ls())
library(qs2)
library(Seurat)
library(SCENIC)
library(dplyr)
library(patchwork)
library(ggplot2) 
library(stringr)
library(plotly)

rss = qs_read("rss.qs2")
colnames(rss) <- paste0("cluster", colnames(rss))

rssPlot <- plotRSS(rss)
rssPlot$plot
plotly::ggplotly(rssPlot$plot)

# 获取每个 cluster 前 10 个 TF 的名字
top0 <- rownames(rss[order(-rss[, "cluster0"]), ])[1:10]
# > top0
#  [1] "NR2F1(+)"  "GSC(+)"    "TRAF4(+)"  "TGIF1(+)"  "IRF7(+)"  
#  [6] "GTF2F1(+)" "CEBPA(+)"  "TCF4(+)"   "BATF3(+)"  "FOXC2(+)" 
top1 <- rownames(rss[order(-rss[, "cluster1"]), ])[1:10]
# > top1
#  [1] "FOXP1(+)"   "ZNF343(+)"  "BHLHE41(+)" "SRF(+)"     "EGR4(+)"   
#  [6] "FOXF2(+)"   "HLF(+)"     "ESRRA(+)"   "MXI1(+)"    "FOXN3(+)" 
top2 <- rownames(rss[order(-rss[, "cluster2"]), ])[1:10]
# > top2
#  [1] "PITX1(+)" "FOXF2(+)" "PHF8(+)"  "ERG(+)"   "CREB5(+)" "HAND2(+)"
#  [7] "PML(+)"   "NR2F6(+)" "HSF1(+)"  "KLF6(+)" 
# 合并去重，得到所有要展示的 TF 名单
top_TFs <- unique(c(top0, top1, top2))
# 提取这些 TF 的 RSS 数据
rss_subset <- rss[top_TFs, ]
# 重新排序行，方便图示（可选）
rss_subset <- rss_subset[order(-apply(rss_subset, 1, max)), ]
# 重新绘图，只展示选中的 TF
rssPlot <- plotRSS(rss_subset)
rssPlot$plot


# 分细胞类型展示regulon特异性分数
p1 = plotRSS_oneSet(rss, 'cluster0')
p2 = plotRSS_oneSet(rss, 'cluster1')
p3 = plotRSS_oneSet(rss, 'cluster2')

p4 = p1 | p2 | p3
p4

ggsave("rank_plot.png", dpi=300)