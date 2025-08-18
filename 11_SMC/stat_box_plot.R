library(tidyverse)
library(Seurat)
library(qs2)
library(future)
library(ggpubr)
plan("multisession", workers =5)
options(future.globals.maxSize = 10*1024^3) # 每个线程分配 10 GB 运行内存
rm(list = ls())
setwd("./11_SMC/")
SMC_cluster = qs_read("./SMC_cluster.qs2")

df = SMC_cluster@meta.data
# Step 1: 添加一个新列 big_group，分为 Normal 和 HF
df <- df %>%
  mutate(big_group = ifelse(group == "DM", "DM", "Others"))

# Step 2: 计算每个样本中 decrease 类型的比例
# 假设每行是一个细胞，因此我们统计每个样本中 decrease 的比例
proportion_df <- df %>%
  group_by(orig.ident, big_group) %>%
  summarise(
    total = n(),
    decrease_count = sum(SMC_type == "decrease"),
    decrease_ratio = decrease_count / total,
    .groups = 'drop'
  )
proportion_df$big_group <- factor(proportion_df$big_group, levels = c("Others", "DM"))
# Step 3: 绘制箱线图
ggplot(proportion_df, aes(x = big_group, y = decrease_ratio, fill = big_group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5) +
  labs(x = "", y = "Proportion of decrease fibro cells") +
  theme_pubr() +
  theme(
    legend.position = "none",
    text = element_text(size = 14)
  )
ggsave("./05_fibro_annotation/fibro_compare_boxplot.png", dpi=300)
