dcp = get_deg_all(exp, group_list, ids, entriz = FALSE, adjust=FALSE, logFC_cutoff = 1, pvalue_cutoff = 0.05)
38 down genes,36 up genes

差异基因数量好少，而且感觉质控有问题，但是之前的图不是这样的，是因为使用了limma的normlize函数吗？
不是的，是因为之前的logFC使用了0.5，现在使用了1，还是改成0.585吧，差异基因控制在200~2000比较合适，否则不方便取交集了

不使用 limma normalize，得到 33 down genes,53 up genes，区别不大，因为这个数据集log之后不同样本直接已经比较对其了

但是质控图明显好看了，样本聚类符合预期
因为聚类热图是使用了差异基因，不是全部基因（会很慢、很大），所以差异的阈值会影响热图，在tinyarray不建议太严格，差异基因不能太少

但是感觉这个数据集质量不高啊，使用所有和部分基因绘图（使用不同的logFC），聚类结果不一致，奇怪
