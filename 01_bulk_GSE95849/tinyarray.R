Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(GEOquery)
library(dplyr)
library(tinyarray)
library(tidyverse)

setwd("./01_bulk_GSE95849/")

# GEO数据下载
# 直接用GSE号即可，默认会通过曾老师的GEO中国镜像下载，超级快，不需要代理
# 下载来的gse是一个长度为3的列表，第1个就是表达矩阵，第2个是样本信息（pdata），第3个是GPL信息
geoID = "GSE95849"
gse <- geo_download(geoID) # 如果报错可添加by_annopbrobe = T从官方途径下载。
# gse = getGEO("GSE95849", destdir=".", AnnotGPL = F, getGPL = F)
exp=gse$exp #从名为 geo 的对象中提取名为exp的组件，并将提取的组件赋值给一个新的变量exp。
exp[1:4,1:4]#自行判断是否需要进行log
gse$exp = log2(gse$exp+1)

boxplot(gse$exp)#检查数据是否异常

exp = gse$exp

soft=getGEO(filename ="GPL22448_family.soft")
gpl = soft@dataTable@table

gene = gpl[,c(1,4)]
###exp是矩阵 merge合并必须是数据框的合并
exp2=as.data.frame(exp)
exp.anno=merge(x=gene,y=exp2,by.x=1,by.y=0)
x=gene$Gene_symbol

# 也要注意 gene symbol有些是以-分割的，要注意，但是不要误伤线粒体基因
# 字符串切割
a1=strsplit(x,split = "|",fixed = T) # 注意///前后有空格，精确匹配
##把分割后的字符串的第一个元素提取出来，合并成为一个新的向量
gene.all = sapply(a1,function(x){x[1]})
exp.anno$Gene_symbol=gene.all
###到此为止  我们讲注释的gene symbol文件与表达矩阵都已经准备好了
###########################################################################################
# 进行探针的转化
exp1=exp.anno
colnames(exp1)[2]="gene.all"
# 整理表达矩阵
# duplicate 'row.names' are not allowed
exp2 = distinct(exp1,gene.all,.keep_all = T) # 这里我们取第一个

# ! missing values in 'row.names' are not allowed
exp3 = exp2[!is.na(exp2$gene.all),]
exp3 = exp3 |> filter(gene.all != "previous version conserved probe")
# 为什么只少了一条呢？因为前面去除重复的时候把NA也只保留了一个
rownames(exp3)=exp3$gene.all

###去除第一列和第二列
exp4=exp3[,-c(1,2)]

saveRDS(exp4, paste0(geoID, "_exp.rds"))