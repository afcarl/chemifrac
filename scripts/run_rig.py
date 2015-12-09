
import numpy as np
import networkx as nx
import pandas as pd
import copy
from skbio import DistanceMatrix, OrdinationResults
from skbio.stats.ordination import pcoa
from fractions import Fraction
from scipy.sparse import coo_matrix
import warnings
from rig import rig

meta_file = '../data/coral_meta.txt'
xls_file = '../data/Coral_ChemiFRAC_test.xlsx'

table = pd.read_excel(xls_file, sheetname=1, index_col=0).T
edges = pd.read_excel(xls_file, sheetname=0)
maxID = max([edges['CLUSTERID1'].max(), edges['CLUSTERID2'].max()]) + 1
spm = coo_matrix((edges['Cosine'].values, 
                  (edges['CLUSTERID1'].values, 
                   edges['CLUSTERID2'].values)),
                shape=(maxID, maxID))
coral_nwk = nx.from_scipy_sparse_matrix(spm)
meta_map = pd.read_table('../data/%s' % meta_file) 

small_table = table

dm = pd.DataFrame(columns=meta_map.index,
                  index=meta_map.index)
for i in range(len(meta_map.index)):
    for j in range(i):
        sampIDs = meta_map['#SampleID'].values
        _x, _y = sampIDs[i], sampIDs[j]
        
        x = small_table.loc[_x, :]
        y = small_table.loc[_y, :]
        dm.loc[_x, _y] = rig(coral_nwk, x, y)

dm.to_csv('../results/rig.txt', sep='\t')
dm = pd.read_csv('../results/rig.txt', index_col=0)
dm = dm.loc[meta_map['#SampleID'].values, meta_map['#SampleID'].values]
dm = dm.fillna(0)
dmpc = pcoa(dm + dm.T)
dmpc.samples.index = meta_map['#SampleID'].values
dmpc.write('../results/rig_pc.txt')
