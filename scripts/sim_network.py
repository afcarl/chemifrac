"""
Simulation of metabolic network
"""
from __future__ import division
from frac import frac, sample_dm
from sim import sample_components, draw_graph, merge

import numpy as np
import networkx as nx
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
from scipy.sparse import coo_matrix
res_dir = '../results/simulation'
xls_file = '../data/Coral_ChemiFRAC_test.xlsx'
meta_file = '../data/coral_meta.txt'
sdm = pd.read_excel(xls_file, sheetname=0)

maxID = max([sdm['CLUSTERID1'].max(),  sdm['CLUSTERID2'].max()])+1
dm = coo_matrix((sdm['Cosine'],
                 (sdm['CLUSTERID1'].values,
                  sdm['CLUSTERID2'].values)),
		shape=(maxID, maxID))
dmG = nx.from_scipy_sparse_matrix(scipy.sparse.dok_matrix(dm))

np.random.seed(0)
sampled_graphs = sample_components(dmG,
                                   multitons=2,
                                   singletons=2)
sim_graph = merge(sampled_graphs)

plt.close('all')
fig = draw_graph(dmG, nx.circular_layout)
fig.savefig('%s/network.png' % res_dir)

# Now to build a bucket table with two samples
dm = nx.adjacency_matrix(sim_graph)
# sdm = pd.DataFrame(scipy.sparse.csgraph.dijkstra(dm),
#                    index=sim_graph.node.keys(),
#                    columns=sim_graph.node.keys())
sdm = scipy.sparse.csgraph.dijkstra(dm)

# The scenario where samples don't have any common metabolites
env1 = pd.Series({119796: 1/5,
                  142883: 0,
                  154390: 0,
                  178554: 0,
                  204966: 0,
                  283976: 1/5,
                  284168: 1/5,
                  285464: 1/5,
                  289481: 1/5})
env2 = pd.Series({119796: 0,
                  142883: 1/4,
                  154390: 1/4,
                  178554: 1/4,
                  204966: 1/4,
                  283976: 0,
                  284168: 0,
                  285464: 0,
                  289481: 0})
sim_table=pd.DataFrame({'samp1' : generate_sample(env1),
                        'samp2' : generate_sample(env2)})
sim_dm1 = sample_dm(sim_table.T, sdm)

# The scenario where the two samples share 1 metabolite
env1 = pd.Series({119796: 0,
                  142883: 0,
                  154390: 0,
                  178554: 0,
                  204966: 0,
                  283976: 1/4,
                  284168: 1/4,
                  285464: 1/4,
                  289481: 1/4})
env2 = pd.Series({119796: 1/5,
                  142883: 1/5,
                  154390: 1/5,
                  178554: 1/5,
                  204966: 1/5,
                  283976: 0,
                  284168: 0,
                  285464: 0,
                  289481: 0})
sim_table=pd.DataFrame({'samp1' : generate_sample(env1),
                        'samp2' : generate_sample(env2)})
sim_dm2 = sample_dm(sim_table.T, dm.todok())

# The scenario where the samples share all of the same metabolites
env1 = pd.Series({119796: 0,
                  142883: 0,
                  154390: 0,
                  178554: 0,
                  204966: 0,
                  283976: 1/4,
                  284168: 1/4,
                  285464: 1/4,
                  289481: 1/4})
env2 = pd.Series({119796: 0,
                  142883: 0,
                  154390: 0,
                  178554: 0,
                  204966: 0,
                  283976: 1/4,
                  284168: 1/4,
                  285464: 1/4,
                  289481: 1/4})
sim_table=pd.DataFrame({'samp1' : generate_sample(env1),
                        'samp2' : generate_sample(env2)})
sim_dm3 = sample_dm(sim_table.T, dm.todok())
