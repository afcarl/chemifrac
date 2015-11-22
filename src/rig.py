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


def rig(G, x, y):
    """ Compute the RIG metric on all components.

    Parameters
    ----------
    G : nx.Graph
       A connected graph of weighted edges
    x : list, str
       Set of nodes present in sample x
    y : list, str
       Set of nodes present in sample y

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
    for comp in nx.connected_components(G):
        nodes = set(comp.nodes())
        cost += rig_component(comp, x & nodes, y & nodes)
    return cost

def rig_component(G, x, y):
    """ Compute the RIG metric on a single component.

    Parameters
    ----------
    G : nx.Graph
       A connected graph of weighted edges
    x : list, str
       Set of nodes present in sample x
    y : list, str
       Set of nodes present in sample y

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
    if len(x)==0 and len(y)==0:
        return 0

    cost = 0
    edges = G.edges(data='weight')
    # This normalizes all of the distances by summing over
    # all of the edge weights. This is doubled, since edges
    # can be counted twice
    weight = sum(list(zip(*edges))[2]) * 2

    # If there is no overlap, arbituarily pick two metabolites
    # and append them to the set.  This to address the issue of
    # measuring distance between unshared components
    if len(x & y) == 0:
        x = x | set(list(y)[0])
        y = y | set(list(x)[0])
        cost += 1

    if x.issubset(y):
        alignedG = _subsetted_graph(G, x, y)
    elif y.issubset(x):
        alignedG = _subsetted_graph(G, y, x)
    else:
        alignedG = _non_subsetted_graph(G, x, y)
    W, _ = nx.network_simplex(alignedG.to_directed())
    cost += W / weight
    return cost

def rig_grow(G, x, y):
    """ Find the minimum weighted path cover on a single component.

    This is a special case of path cover finding the minimum weighted path
    cover linking sample x and sample y.  Ultimately, this boils down to the
    mininum flow problem, which can be solved via network simplex algorithm.

    Parameters
    ----------
    G : nx.Graph
       A connected graph of weighted edges
    x : list, str
       Set of nodes present in sample x
    y : list, str
       Set of nodes present in sample y

    Returns
    -------
    float :
       Distance between sample x and sample y

    Note
    ----
    If x or y is None, then 1 will be added to the total distance.
    If they are both None, then the distance will be zero.

    Also, this module as it stands is only dealing with presence/absence

    """
    if x.issubset(y):
        alignedG = _subsetted_graph(G, x, y)
    elif y.issubset(x):
        alignedG = _subsetted_graph(G, y, x)
    else:
        alignedG = _non_subsetted_graph(G, x, y)
    cost, edges = nx.network_simplex(alignedG.to_directed())
    return cost

def _subsetted_graph(G, x, y):
    """ Finds mininum flow cover in the scenario where x is a subset of y

    This is a special case of path cover finding the minimum weighted path
    cover linking sample x and sample y.  Ultimately, this boils down to the
    mininum flow problem, which can be solved via network simplex algorithm.

    Parameters
    ----------
    G : nx.Graph
       A connected graph of weighted edges
    x : list, str
       Set of nodes present in sample x.  x is a subset of y
    y : list, str
       Set of nodes present in sample y.  y is a superset of x

    Returns
    -------
    nx.Graph :
       Aligned graph

    """
    xtotal = len(x)

    inX, outX = x, y - x
    newG = copy.deepcopy(G)
    for node in G.nodes():
        if node in inX:
            newG.node[node]['demand'] = Fraction(1, xtotal)
        else:
            newG.node[node]['demand'] = Fraction(-1, xtotal)
    return newG

def _non_subsetted_graph(G, x, y):
    """ Finds mininum flow cover in the scenario where x is not a subset of y

    This is a special case of path cover finding the minimum weighted path
    cover linking sample x and sample y.  Ultimately, this boils down to the
    mininum flow problem, which can be solved via network simplex algorithm.

    Parameters
    ----------
    G : nx.Graph
       A connected graph of weighted edges
    x : list, str
       Set of nodes present in sample x.  x is a subset of y
    y : list, str
       Set of nodes present in sample y.  y is a superset of x

    Returns
    -------
    nx.Graph :
       Aligned graph

    """
    Gx = copy.deepcopy(G)
    Gy = copy.deepcopy(G)
    num_shared = len(x & y)
    num_uniquex = len(x - y)
    num_uniquey = len(y - x)
    xtotal = len(x)
    ytotal = len(y)
    for node in G.nodes():
        if node in x and node not in y:
            Gx.node[node]['demand'] = Fraction(1, xtotal)
            Gy.node[node]['demand'] = Fraction(-1, num_uniquex)
        elif node not in x and node in y:
            Gx.node[node]['demand'] = Fraction(-1, num_uniquey)
            Gy.node[node]['demand'] = Fraction(1, ytotal)
        elif node in x and node in y:
            Gx.node[node]['demand'] = Fraction(1, xtotal)
            Gy.node[node]['demand'] = Fraction(1, ytotal)
        else:
            Gx.node[node]['demand'] = 0
            Gy.node[node]['demand'] = 0
    alignedG = nx.union(Gx, Gy, rename=('x-','y-') )
    return alignedG
