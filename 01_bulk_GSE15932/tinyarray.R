rm(list = ls())
library(tinyarray)
library(tidyverse)

getwd()
setwd("./DM_ED/01_bulk_GSE15932/")

# GEO数据下载
# 直接用GSE号即可，默认会通过曾老师的GEO中国镜像下载，超级快，不需要代理
# 下载来的gse是一个长度为3的列表，第1个就是表达矩阵，第2个是样本信息（pdata），第3个是GPL信息
geoID = "GSE15932"
gse <- geo_download(geoID) # 如果报错可添加by_annopbrobe = T从官方途径下载。

exp=gse$exp #从名为 geo 的对象中提取名为exp的组件，并将提取的组件赋值给一个新的变量exp。
exp[1:4,1:4]#自行判断是否需要进行log
gse$exp = log2(gse$exp+1)
exp = exp[,c(9:16,25:32)]

boxplot(gse$exp)#检查数据是否异常
# 数据分组
pd = gse$pd
pd = pd[c(9:16,25:32),]
glimpse(pd)
Group=ifelse(str_detect(pd$title,"diabetes"), # 自行寻找分组信息所在的组合分组的情况
             "diabetes",
             "Control")
Group = factor(Group,levels = c("Control", "diabetes"))
Group
# Group = data.frame(
#   sample = colnames(exp),
#   Group = Group
# )
# write.csv(Group, "Group.csv",row.names = FALSE)

# ID转换
ids <- AnnoProbe::idmap(gse$gpl) # 配合AnnoProbe

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
exp1 = log2(exp+1)
saveRDS(exp1, paste0(geoID, "_exp.rds"))
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

# 差异分析及可视化
exp0 = gse$exp
exp0 = exp0[,c(9:16,25:32)]
dcp = get_deg_all(exp0,Group,ids, entriz = FALSE, adjust=FALSE, logFC_cutoff = 0.5) # 根据提供表达矩阵、分组信息和探针注释，返回差异分析结果
table(dcp$deg$change)
write.csv(dcp$deg, paste0(geoID, "_deg.csv"), row.names=FALSE)
dcp$plots
ggsave(paste0(geoID, "-degplots.png"),width = 15,height = 5, bg="white", dpi=300)
