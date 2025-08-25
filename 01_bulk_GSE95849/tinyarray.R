rm(list = ls())
library(tidyverse)
library(patchwork)
library(GEOquery)
library(tinyarray)
library(limma)

getwd()
setwd("./01_bulk_GSE95849/")

# GEO array 数据下载
# 直接用GSE号即可，默认会通过曾老师的GEO中国镜像下载，超级快，不需要代理
geoID = "GSE95849"
file_path <- paste0("./", geoID, "_eSet.Rdata")

# 检查文件是否存在
if (file.exists(file_path)) {
  # 文件存在，直接加载
  load(file_path)
  message("文件已存在，直接加载数据")
} else {
  # 文件不存在，进行下载
  gse <- geo_download(geoID) # 如果报错可添加by_annopbrobe = T从官方途径下载
  message("文件下载完成并已保存")
}


# 数据预处理


# 表达矩阵
gse = gset[[1]]
exp = exprs(gse) # matrix, 从名为 geo 的对象中提取名为exp的组件，并将提取的组件赋值给一个新的变量exp。
boxplot(exp) # 用boxplot来看数据的分布非常重要

exp[1:4,1:4] # 自行判断是否需要进行log
exp = log2(exp + 1)
boxplot(exp) # 检查数据分布

exp = normalizeBetweenArrays(exp)
boxplot(exp) # 检查数据分布

exp = exp[,1:12]
saveRDS(exp, paste0(geoID, "_normalized_exp.rds"))

# 表型信息
pd = pData(gse)
glimpse(pd)
pd = pd[1:12,]
group_list = ifelse(str_detect(pd$title,"DM"), # 自行寻找分组信息所在的组合分组的情况
             "DM",
             "Control")
group_list = factor(group_list,levels = c("Control", "DM")) # control为参考组，在前
group_list
# Group = data.frame(
#   sample = colnames(exp),
#   Group = group_list
# )
# write.csv(Group, "Group.csv",row.names = FALSE)



# ID准备
# 生成文件路径
ids_file <- paste0(gse@annotation, "_ids.rds")

# 检查文件是否存在
if (file.exists(ids_file)) {
  # 文件存在，直接读取
  ids <- readRDS(ids_file)
  message("ID映射文件已存在，直接加载")
} else {
  # 文件不存在，进行生成并保存
  ids <- AnnoProbe::idmap(gse@annotation)
  saveRDS(ids, file = ids_file)
  message("ID映射文件生成完成并已保存")
}


## 自己制作id文件
# 如果不能使用AnnoProbe找到，则自己制作ids
# ids为一个dataframe，第一列为probe_id，第二列为symbol

### 读取平台文件
soft = getGEO(filename ="GPL22448_family.soft")
gpl = soft@dataTable@table
gene = gpl[,c(1,4)]

exp2=as.data.frame(exp) # exp是矩阵 merge合并必须是数据框的合并
exp.anno=merge(x=gene,y=exp2,by.x=1,by.y=0)
x=gene$Gene_symbol

# 字符串切割, 也要注意 gene symbol有些是以-分割的
a1=strsplit(x,split = "|",fixed = T) # 注意///前后有空格，精确匹配
## 把分割后的字符串的第一个元素提取出来，合并成为一个新的向量
gene.all = sapply(a1,function(x){x[1]})
exp.anno$Gene_symbol=gene.all

ids = exp.anno[c(1, 2)]
colnames(ids) = c("probe_id", "symbol")
ids = distinct(ids, symbol, .keep_all = T)
ids = ids[!is.na(ids$symbol),]
ids = ids |> filter(symbol != "previous version conserved probe")
saveRDS(ids, paste0(gse@annotation, "_ids"))
# ID转换
# 因为探针和基因不是一一对应的关系，会存在多个探针对应同一个基因的情况
# 首先，多个探针对应同一个基因的问题是必须要处理的，因为最终的分析是以基因为单位，一个基因在一个样本里只能有一个表达量，在一次差异分析中只能有一个logFC，不能有多个啊。
# 全部的处理的方法有三种：
# 1.随机去重，保留任何一个探针都可以
# 2.保留行和/行平均值最大的探针
# 3.取多个探针的平均值
# 采取这三种方法都是可以的，没有标准答案

# 1.随机去重
# 因为这些探针也没什么主次和顺序，所以直接用R语言里的去重复函数搞一下就可以。trans_array这个函数采用的就是直接去重，这个例子里就会保留第一个探针。
exp1 <- trans_array(exp, ids)
saveRDS(exp1, paste0(geoID, "_anno_exp.rds"))
# 2.保留行和/行平均值最大的探针
# ids = ids[ids$probe_id %in% rownames(exp),]
# exp = exp[ids$probe_id,]
# identical(ids$probe_id,rownames(exp)) #他俩的顺序已经相同了
# ids = ids[order(rowSums(exp),decreasing = T),] #ids可以使用rowSums求和的顺序来排序，从大到小排，最大的在上面
# exp2 = trans_array(exp,ids)

# # 3.求平均值
# exp3 = exp[ids$probe_id,]
# rownames(exp3) = ids$symbol
# exp3[1:4,1:4]
# exp3 = limma::avereps(exp3)
# saveRDS(exp3, file = "exp3.rds")


# 差异分析及可视化，一步完成
dcp = get_deg_all(exp, group_list, ids, entriz = FALSE, adjust=FALSE, logFC_cutoff = 0.585, pvalue_cutoff = 0.05) # 根据提供表达矩阵、分组信息和探针注释，返回差异分析结果
write.csv(dcp$deg, paste0(geoID, "_deg.csv"), row.names=FALSE)
dcp$plots # 注意质量控制 
ggsave(paste0(geoID, "_plots.png"),width = 15,height = 5, bg="white", dpi=300)

# 左边的热图，说明我们实验的两个分组，normal和npc的很多基因表达量是有明显差异的
# 中间的PCA图，说明我们的normal和npc两个分组非常明显的差异
# 右边的火山图展示了差异分析的结果

# 差异分析及可视化，分步(方便自定义调整)
pvalue_cutoff = 0.05
logFC_cutoff = 1
adjust = FALSE
entriz = FALSE # 函数返回的结果数据框 deg 中是否会多出一列 ENTREZID。
species = "human"
symmetry = TRUE
my_genes = NULL # genes for pheatmap
lab = NA
show_rownames = FALSE
cluster_cols = TRUE
n_cutoff = 2
annotation_legend = FALSE

deg <-  get_deg(exp,group_list,ids,
                logFC_cutoff=logFC_cutoff,
                pvalue_cutoff=pvalue_cutoff,
                adjust = adjust,
                entriz = entriz,
                species = species)
cgs = get_cgs(deg)

volcano_plot = draw_volcano(deg,pkg=4,
                            lab =lab,
                            pvalue_cutoff = pvalue_cutoff,
                            logFC_cutoff=logFC_cutoff,
                            adjust = adjust,
                            symmetry = symmetry)

pca_plot = draw_pca(exp,group_list)
heatmap = draw_heatmap2(exp,group_list,deg,my_genes,
                        show_rownames = show_rownames,
                        n_cutoff = n_cutoff,
                        cluster_cols = cluster_cols,
                        annotation_legend=annotation_legend)
if(as.numeric(grDevices::dev.cur())!=1) grDevices::graphics.off()

wrap_plots(heatmap,pca_plot,volcano_plot)+plot_layout(guides = 'collect')
ggsave(paste0(geoID, "_plots.png"),width = 15,height = 5, bg="white", dpi=300)