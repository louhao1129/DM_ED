import os,sys
import loompy as lp
import numpy as np
import pandas as pd
os.getcwd()
os.chdir("./14_pyscenic\\")

counts=pd.read_csv("./result/counts.csv",index_col = 0, header = 0)
features = counts.index.tolist()
barcodes = counts.columns.tolist()
#loom行属性：用于存储与基因相关的元数据。这些属性可以包括基因名称、基因长度、基因所属的染色体等信息。
#loom列属性：存储与细胞相关的元数据。例如细胞类型等信息。
row_attrs = {'Gene':np.array(features)}
col_attrs = {'CellID':np.array(barcodes)}
lp.create("./result/counts.loom", counts.values,row_attrs,col_attrs)
