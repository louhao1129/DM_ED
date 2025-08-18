-   GSE206528

    -   来自三名勃起正常的男性和五名器质性勃起功能障碍 （ED） 患者的海绵体单细胞转录组。ED包括 2 名糖尿病勃起功能障碍 （DMED） 患者和 3 名非糖尿病性勃起功能障碍 （non-DM） 患者
    -   根据以下阈值参数进一步过滤细胞：表达基因的总数，800-7000;总 UMI 计数，0-30,000;和表达的线粒体基因比例，\<10%。使用 Seurat 包装中的 IntegrateData 函数进行批量校正。
    -   合并的 Seurat 对象通过主成分分析 （PCA） 进行缩放和分析。前 20 台 PC 还用于进行聚类和执行 t 分布随机邻域嵌入 （tSNE） 降维。使用 Seurat 包中的 FindClusters 函数（分辨率参数设置为 0.5）对单元格进行聚类。为了进一步分析每个集群，我们将它们隔离并再次执行上述两个步骤以获取子集群信息。
    -   对于每个主要簇，通过 PECAM1 和 VWF 的高表达来识别 EC 簇;通过 LUM、COL1A1 和 PDGFRA 的表达鉴定 FB 簇;通过 RGS5 和 KCNJ8 高表达，MYOCD 低表达来鉴定 PC 簇;通过高表达 ACTA2 和 MYH11 鉴定 SMC 簇;通过 S100B 和 MPZ 的高表达来鉴定 SWC 簇;通过 CD163 和 CD68 的高表达来鉴定 MAC 簇;通过 CD3D 和 CD3E 的高表达鉴定 T 簇。
    -   使用 Seurat 函数 FindAllMarkers（test.use = wilcox;min.pct = 0.1;logfc.threshold = 0.25）根据归一化 UMI 计数鉴定差异表达基因（DEG）。除非另有说明，否则每个选定子聚类中的 DEG 是根据该子聚类与数据集其余部分之间的比较计算的。使用“WebGestalt \[网站 <http://www.webgestalt.org%5D”进行> GO 分析，选择过表分析（ORA）或基因集富集分析（GSEA）作为感兴趣的方法，在功能数据库中仅选择生物过程。使用基于 DEG 的 log2 （FC） 和 P 值的 Ingenuity Pathway Analysis （IPA） 软件进行通路分析。本研究中使用的其他基因集列在补充数据 15 中。