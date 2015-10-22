"""
The RIG metric

Retract
Insert
Grow
"""
import numpy as np
import networkx as nx
import pandas as pd
from scipy.sparse.csgraph import dijkstra


def growth(G):
    """
    Performs the growth step on a single component

    Parameters
    ----------
    G : nx.Graph
       Connected graph.  Each edge has a weight and each node has a
       pair of values - each value corresponding to presence/absence in
       each sample

    Returns
    -------
    modG : nx.Graph
       Modified component
    d : float
       distance between the two samples calculated on the component
    """
    shortest_paths = pd.DataFrame(dijkstra(nx.to_scipy_sparse_matrix(G)),
                                  index=G.node.keys(),
                                  columns=G.node.keys())
    # Iterate through all nodes in sample 1
    # For each node in sample 1, find the shortest path to any node in sample 2

    # Repeat this for sample 2

def _grow_sample(G, dm, idx=0):
    """
    Expands sample [idx] to all of the metabolites contained in sample[!idx]

    Parameters
    ----------
    G : nx.Graph
       Connected graph.  Each edge has a weight and each node has a
       pair of values - each value corresponding to presence/absence in
       each sample
    dm: pd.DataFrame
       Shortest path distance matrix
    idx: bool
       Index of sample to modify
    """
    d, dest = _shortest_distance(G, dm, source)



def _shortest_distance(G, dm, source, idx):
    """
    Finds the shortest path from source node
    to any node contained in the other sample

    Parameters
    ----------
    G : nx.Graph
       Connected graph.  Each edge has a weight and each node has a
       pair of values - each value corresponding to presence/absence in
       each sample
    dm: pd.DataFrame
       Shortest path distance matrix
    source: str
       Source node to find path connecting to nodes contained in the
       other sample
    idx: bool
       Index of sample to modify

    Returns
    -------
    d : float
       Distance to nearest source contained in the other sample
    dest : str
       Name of the closest destination node found in the other sample
    """
    assert len(G.node[source]['presence'])==2
    other_idx = np.logical_not(idx)
    other_keys = [n for n in G.node.keys()
                  if G.node[n]['presence'][other_idx]==1]
    paths = pd.Series([dm.loc[source, i]
                       for i in other_keys],
                      index=other_keys)
    d, dest = paths.min(), paths.argmin()
    return d, dest
