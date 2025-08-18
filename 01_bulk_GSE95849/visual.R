library(tidyverse)                         
library(ggrepel)
library(ggsci)  
library(ggsignif)
library(ggpubr)

data = readRDS("./GSE95849_deg.rds")
data$significant="stable"
log2(1.5)
data$significant[data$logFC>=0.585 & data$P.Value <0.05]="up" # 1.5倍，一般差异基因在500~2000之间
data$significant[data$logFC<= -0.585 & data$P.Value <0.05]="down"

ggplot(data,aes(logFC,-1*log10(P.Value)))+xlim(-2,2)+ylim(0,6)+ # x、y轴范围要视情况调整
  geom_point(aes(color=significant),size=0.8)+theme_classic()+
  scale_color_manual(values = c("#2a9d8f","#8d718cff","#f8961e"))+
  geom_hline(yintercept = 1.3,linetype=4,size=0.3)+
  geom_vline(xintercept = c(-0.585,0.585),linetype=4,size=0.3)+
  theme(title=element_text(size = 18),text = element_text(size=18))+
  labs(x="log2(foldchange)",y="-log10(p_value)")

select.FPKM <- data$AveExpr > 5 
table(select.FPKM)

select.log2FC <- abs(data$logFC) >0.5
table(select.log2FC)
select.qval <- (data$adj.P.Val< 0.05)
table(select.qval)

select.vec=(select.FPKM & select.log2FC & select.qval) 
table(select.vec)

degs.list=as.character(rownames(data))[select.vec]

label.deg=sample(degs.list,20) 
# 这里是随机抽样，也可以选取自己感兴趣的基因
# label.deg = c()
p=ggplot(data,aes(logFC,-1*log10(P.Value)))+xlim(-2,2)+ylim(0,6)+
  geom_point(aes(color=significant),size=0.8)+theme_classic()+
  scale_color_manual(values = c("#2a9d8f","#EE7AE9","#f8961e"))+
  geom_hline(yintercept = 1.3,linetype=4,size=0.3)+
  geom_vline(xintercept = c(-0.5,0.5),linetype=4,size=0.3)+
  theme(title=element_text(size = 18),text = element_text(size=18))+
  labs(x="log2(foldchange)",y="-log10(p_value)")
data_selected <- data[label.deg,]
p + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))
###########################
data_selected <- data["PCNX2",]
p + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))


####################################################3
##画个热图

annotation_col1 = data.frame(
  Database =c(rep("GEO",12)),
  CellType =c(rep("control",6),rep("treatment",6)) 
)
rownames(annotation_col1)=colnames(exp.5)
exp.6=as.data.frame(exp.5)
exprSet.map=exp.6[label.deg,]

pheatmap::pheatmap(exprSet.map, #热图的数据
                   cluster_rows = F,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   annotation_col =annotation_col1,
                   show_colnames=F,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("blue", "white","red"))(100))
####聚类
pheatmap::pheatmap(exprSet.map, #热图的数据
                   cluster_rows =T,#行聚类
                   cluster_cols =F,#列聚类，可以看出样本之间的区分度
                   annotation_col =annotation_col1,
                   show_colnames=F,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("blue", "white","red"))(100))

######################################################
###箱线图
label.deg[1:10]
exp.x=exp.6["TTC32",]

z=t(exp.x)
z=as.data.frame(z)
z$type=c(rep("control",6),rep("treatment",6))

colnames(z)
library(ggpubr)
library(ggsignif)
library(ggplot2)


ggboxplot(z,x="type",y="TTC32",
          width = 0.6,fill="type",
          notch = T,palette = c("#00AFBB", "red","#E7B800"),
          add = "jitter",shape="type")

###加个p值
p=ggboxplot(z,x="type",y="TTC32",
          width = 0.6,fill="type",
          notch = T,palette = c("#00AFBB", "red","#E7B800"),
          add = "jitter",shape="type")

p + stat_compare_means(aes(group =type))

####################################################
###箱线图
ggplot(z,aes(type,TTC32,fill=type)) +
  geom_boxplot()



#######匹配一下，注意pdata看哪些样本来自一个个体， 配对箱线图
z$pairinfo=pairinfo=rep(1:6,2)
ggplot(z, aes(type,TTC32,fill=type)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5) +
  geom_line(aes(group=pairinfo), colour="black", linetype="11") +
  xlab("") +
  ylab(paste("Expression of ","TTC32"))+
  theme_classic()+
  theme(legend.position = "none")

##########################################################################################
#####小提琴

ggplot(z,aes(type,TTC32,fill=type)) +
  geom_violin()+ geom_jitter(shape=16, position=position_jitter(0.2))


compaired <- list(c("control", "treatment"))



colnames(z)
ggplot(z,aes(type,TTC32,fill=type)) +
  geom_violin()+
  geom_signif(comparisons =compaired ,step_increase = 0.1,map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),test =t.test,size=2,textsize = 6)+
  geom_jitter(shape=16, position=position_jitter(0.2))

help(package="ggsci")


