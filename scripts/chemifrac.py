from __future__ import division
from frac import frac, sample_dm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re
from skbio import TreeNode
from StringIO import StringIO
from biom import load_table
import pandas as pd
import networkx as nx
import scipy
import copy
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.distance import DistanceMatrix
from skbio.diversity.beta import pw_distances
from skbio.stats.ordination import pcoa
from scipy.sparse.csgraph import dijkstra
from scipy.spatial.distance import pdist

# Calculates chemifrac distance
xls_file = '../data/Coral_ChemiFRAC_test.xlsx'
meta_file = '../data/coral_meta.txt'
table = pd.read_excel(xls_file, sheetname=1, index_col=0)
meta_map = pd.read_table(meta_file)

p = meta_map.loc[meta_map['Organism'] == 'Porites', :]
b = meta_map.loc[meta_map['Organism'] == 'Black Nasty', :]
i = meta_map.loc[meta_map['Organism'] == 'Interaction', :]
meta_map = pd.concat([p, b, i])

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

# full_dm = copy.deepcopy(dm)
# full_dm = full_dm.fillna(0)
# dmG = nx.from_scipy_sparse_matrix(scipy.sparse.dok_matrix(full_dm))
# connected_comps = list(nx.connected_components(dmG))
# connected_graphs = list(nx.connected_component_subgraphs(dmG))
# for i, g in enumerate(connected_graphs):
#     nx.write_gml(g, '../results/graphs/%d.gml'%i)


# dm = dm.fillna(1)
# dm = dm.iloc[:500, :500]
otu_table = table.T
otu_table = otu_table.loc[:, dm.index]
dm = dm.fillna(np.inf)
cosine = dm

cosine = cosine.reindex_axis(sorted(otu_table.columns), axis=0)
cosine = cosine.reindex_axis(sorted(otu_table.columns), axis=1)

shortest = dijkstra(dm.values)
shortest = pd.DataFrame(shortest,
                        columns=dm.index, index=dm.columns)
shortest = shortest.reindex_axis(sorted(otu_table.columns), axis=0)
shortest = shortest.reindex_axis(sorted(otu_table.columns), axis=1)
# otu_table = table.T
otu_table = otu_table.reindex_axis(sorted(otu_table.columns), axis=1)

# Uses an idea similar to simrank
graph_dm = (otu_table>0).dot(cosine).dot((otu_table>0).T)
graph_dm.to_csv('../results/simrank.txt', '\t')
# Uses Aitchison distance
# samples = ['CF31_A', u'CF31_B', u'CF141_A', u'CF141_B', u'Tuni', u'Bry']
dm = cosine.values
dm[dm==np.inf]=0
mat = otu_table.values
mat = multiplicative_replacement(mat)
graph_dm = connected_dm(mat, dm)
graph_dm += graph_dm.T
samples = otu_table.index
graph_dm = pd.DataFrame(graph_dm,
                        index=samples,
                        columns=samples)
graph_dm.to_csv('../results/aitchison.txt', '\t')

# Read in graph_dm
graph_dm = pd.read_csv('../results/unconnected_aitchison.txt',
                       sep='\t', index_col=0)
# table = pd.read_table('../data/skinmap_chemiFrac_test.txt',
#                        sep='\t', index_col=0)
graph_dm.index = table.columns
graph_dm.columns = table.columns
# _dm = pw_distances('braycurtis', table.values, table.index.values)
# _dm.write('../results/braycurtis.txt')
_dm = DistanceMatrix(graph_dm.values + graph_dm.values.T)
_dm.ids = graph_dm.index
pcoa_v = pcoa(_dm)

fig = plt.figure(3)
plt.plot(pcoa_v.samples['PC1'],
         pcoa_v.samples['PC2'], 'ob')
# plt.plot(pcoa_v.eigvecs[not_stressed, 0],
#          pcoa_v.eigvecs[not_stressed, 1],
#          'o', c='#FFFFFF', label='Before stress')
# plt.plot(pcoa_v.eigvecs[stressed, 0],
#          pcoa_v.eigvecs[stressed, 1],
#          'o', c='#999999', label='After stress')
# plt.legend(loc=3)
#plt.title('Weighted Aitchison on Coral data')
#fig.savefig('../results/coral_chemifrac.png')

pcoa_v.write('../results/coral_unconnected_aitchison_pcoa.txt')
