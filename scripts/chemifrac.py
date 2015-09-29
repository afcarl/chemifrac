from __future__ import division
from frac import frac, sample_dm, connected_dm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re
from skbio import TreeNode
from StringIO import StringIO
from biom import load_table
import pandas as pd

from skbio.stats.composition import multiplicative_replacement
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import pcoa
from scipy.sparse.csgraph import dijkstra

# Calculates chemifrac distance
xls_file = '../data/Coral_ChemiFRAC_test.xlsx'
table = pd.read_excel(xls_file, sheetname=1, index_col=0)
sdm = pd.read_excel(xls_file, sheetname=0)
ids = list(set(sdm['CLUSTERID1']) | set(sdm['CLUSTERID2']))
dm = pd.DataFrame(index=ids, columns=ids)
for i in range(len(sdm)):
    x, y, d = sdm.ix[i, :]
    x, y = int(x), int(y)
    dm.loc[x, y] = 1-d
    dm.loc[y, x] = 1-d
    dm.loc[x, x] = 0
    dm.loc[y, y] = 0

dm = dm.fillna(np.inf)
# cosine = dm
otu_table = table.T
# cosine = cosine.reindex_axis(sorted(otu_table.columns), axis=0)
# cosine = cosine.reindex_axis(sorted(otu_table.columns), axis=1)

shortest = dijkstra(dm.values)
shortest = pd.DataFrame(shortest,
                        columns=ids, index=ids)
shortest = shortest.reindex_axis(sorted(otu_table.columns), axis=0)
shortest = shortest.reindex_axis(sorted(otu_table.columns), axis=1)

otu_table = otu_table.reindex_axis(sorted(otu_table.columns), axis=1)


# Uses an idea similar to simrank
# graph_dm = (otu_table>0).dot(cosine).dot((otu_table>0).T)
# graph_dm.to_csv('../results/simrank.txt', '\t')
# Uses Aitchison distance
# samples = ['CF31_A', u'CF31_B', u'CF141_A', u'CF141_B', u'Tuni', u'Bry']

samples = otu_table.index
mat = otu_table.values
mat = multiplicative_replacement(mat)
# dmG = nx.Graph(shortest.values)
dm = shortest.values
dm[dm == np.inf] = 0
dm[np.isnan(dm)] = 0
# dmG = nx.Graph(scipy.sparse.dok_matrix(dm))
# dmG = nx.from_scipy_sparse_matrix(scipy.sparse.dok_matrix(dm))

graph_dm = connected_dm(mat, shortest.values)
graph_dm += graph_dm.T
graph_dm = pd.DataFrame(graph_dm,
                        index=samples,
                        columns=samples)
graph_dm.to_csv('../results/unconnected_aitchison.txt', '\t')
