"""
The RIG metric

Retract
Insert
Grow
"""
import numpy as np
import networkx as nx
import pandas as pd
import copy
from scipy.sparse.csgraph import dijkstra
from fractions import Fraction
from sklearn.metrics.pairwise import pairwise_distances
from functools import partial

def rig_pairwise(G, X):
    rig_func = partial(rig, G)
    # Note cannot work with pandas
    return pairwise_distances(X, metric=rig_func)


def rig(G, x, y, res=1e-9):
    """ Compute the RIG metric on all components.

    Parameters
    ----------
    G : nx.Graph
       A connected graph of weighted edges
    x : pd.Series
       Vector of nodes and their abundance in sample x
    y : pd.Series
       Vector of nodes and their abundance in sample y

    Returns
    -------
    float :
       Distance between sample x and sample y

    Note
    ----
    If x or y is None, then 1 will be added to the total distance.
    If they are both None, then the distance will be zero.

    """
    cost = 0
    _G = copy.deepcopy(G)
    # This converts all of the weights to integers
    for u,v,d in _G.edges(data=True):
        d['weight'] = int(d['weight'] / res)

    # This calculates the largest edge set to offset the insertion cost.
    weights = []
    for comp in nx.connected_component_subgraphs(_G):
        edges = comp.edges(data='weight')
        if len(edges) > 0:
            weights.append(sum(list(zip(*edges))[2]))
    maxW = max(weights) + 1

    for comp in nx.connected_component_subgraphs(_G):
        nodes = set(comp.nodes())
        subx = x[nodes & set(x.keys())]
        suby = y[nodes & set(y.keys())]

        # Perform closure so that abundances within
        # each component add up to one
        if len(subx) > 0:
            subx = subx / subx.sum()
        if len(suby) > 0:
            suby = suby / suby.sum()
        c = rig_component(comp, subx, suby, maxW)
        cost += c
    return (cost)*res

def rig_component(G, x, y, maxW):
    """ Compute the RIG metric on a single component.

    Parameters
    ----------
    G : nx.Graph
       A connected graph of weighted edges
    x : pd.Series
       Vector of nodes and their abundance in sample x
       within the connected graph G
    y : pd.Series
       Vector of nodes and their abundance in sample y
       within the connected graph G
    maxW : float
       The cost of an insertion step

    Returns
    -------
    float :
       Distance between sample x and sample y

    Note
    ----
    If x or y is None, then 1 will be added to the total distance.
    If they are both None, then the distance will be zero.
    Also, the weights of the edges must be contained in `'weight'`.
    """

    # Both samples don't contain any metabolites
    if len(x)==0 and len(y)==0:
        return 0

    # Both samples contains the exact same metabolites
    # Note that this will change once we start adding weights
    # due to abundance.
    def equal(x, y):
        if len(x) != len(y):
            return False
        x = x.sort_values(inplace=False)
        y = y.sort_values(inplace=False)

        d = np.allclose(x.values.astype(np.float),
                        y.values.astype(np.float))
        if  d and np.all(x.index==y.index):
            return True
        return False

    if equal(x, y):
        return 0

    # The component being analyzed has only 1 node.
    # So the networkx simplex algorithm cannot be run
    # since there are no edges.
    if len(G.nodes()) == 1:
        return maxW

    cost = 0
    edges = G.edges(data='weight')

    # If there one of the samples doesn't have any metabolites
    # on the component, arbituarily pick two metabolites
    # and append them to the set.  This to address the issue of
    # measuring distance between unshared components.
    if len(x) == 0:
        x = pd.Series({y.index[0]: Fraction(1, 1)})
        weight = maxW
    elif len(y) == 0:
        y = pd.Series({x.index[0]: Fraction(1, 1)})
        weight = maxW
    else:
        weight = 0
    _G = copy.deepcopy(G)

    xarr = pd.Series({n:(x[n] if n in x else 0) for n in G.nodes()})
    yarr = pd.Series({n:(y[n] if n in y else 0) for n in G.nodes()})

    for node in _G.nodes():
        _G.node[node]['demand'] = xarr[node] - yarr[node]

    W, _ = nx.network_simplex(_G.to_directed())

    cost += W + weight
    return cost
