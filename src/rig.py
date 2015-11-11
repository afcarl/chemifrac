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
    If they are both None, then the distance will be zero

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
    cost += rig_grow(G, x, y) / weight
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
    # Create a pair of aligned graphs
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
            Gx.node[node]['demand'] = Fraction(-1,  num_uniquey)
            Gy.node[node]['demand'] = Fraction(1, ytotal)
        elif node in x and node in y:
            Gx.node[node]['demand'] = Fraction(1, xtotal)
            Gy.node[node]['demand'] = Fraction(1, ytotal)
        else:
            Gx.node[node]['demand'] = 0
            Gy.node[node]['demand'] = 0
    alignedG = nx.union(Gx, Gy, rename=('x-','y-') )

    cost, edges = nx.network_simplex(alignedG.to_directed())
    return cost

def insertion(Gs):
    """
    Performs the insertion step on a whole bunch of components

    Parameters
    ----------
    Gs : iterable of nx.Graph
       Connected components.  Each edge has a weight and each node has a
       pair of values - each value corresponding to presence/absence in
       each sample

    Returns
    -------
    modG : nx.Graph
       Modified component
    d : float
       distance between the two samples calculated on the component
    """
    pass

# TODO
# This will need to be recursive
# And will also have to act globally (careful about the greedy approach)
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
    d = 0
    modG = copy.deepcopy(G)
    for idx in [0, 1]:
        for source in G.node.keys():
            # Check both samples simulataneously

            # If the source node isn't present, move on
            if G.node[source]['presence'][idx] == 0: continue
            _d, dest = _shortest_distance(G, dm, source, idx=idx)
            if _d == None: continue
            # If the destination node is present, move one
            if modG.node[dest]['presence'][idx] == 1: continue
            # Set presence to True, hence grow
            modG.node[dest]['presence'][idx] = 1
            d += _d
    return modG, d


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

    Note
    ----
    Returns None, None if no there aren't any other nodes
    present in the other sample that isn't included in the current
    sample
    """
    assert len(G.node[source]['presence'])==2
    other_idx = np.logical_not(idx)
    other_keys = [n for n in G.node.keys()
                  if (G.node[n]['presence'][other_idx]==1
                      and n!=source)]
    paths = pd.Series([dm.loc[source, i]
                       for i in other_keys],
                      index=other_keys)
    if len(paths)==0:
        return None, None

    d, dest = paths.min(), paths.argmin()
    return d, dest

def min_weighted_path_cover(G, x, y):
    """ Find the minimum weighted path cover

    This is a special case of path cover finding the minimum weighted path cover
    linking sample x and sample y

    Parameters
    ----------
    G : nx.Graph
       A connected graph of weighted edges
    x : list, str
       List of nodes present in sample x
    y : list, str
       List of nodes present in sample y

    Returns
    -------
    E : set, nx.edge
       Set of edges within the mininum edge cover
    weight : int
       Total weight of the edges in the edge cover
    """
    E = set() # The edges contained in the optimal edge cover
    xnodes = set(x)
    ynodes = set(y)

    #
    edges = G.edges(data='weight')
    new_edges = list()

    # Filter out all edges not connecting x and y
    edges = _filter_edges(edges, xnodes, ynodes)

    # Sort all of the edges by weight
    edges = sorted(new_edges, key=lambda x: x[2])

    # Since the minimum edge is guaranteed to be in the optimal solution
    # add it to the solution edge set
    minE = edges.pop(0)
    E.add(minE)

    # The remaining nodes that have not been observed yet in the edge set
    Rx, Ry = copy.deepcopy(xnodes), copy.deepcopy(ynodes)

    while len(Rx) > 0 and len(Ry) > 0:
        # Three possibilities
        # 1. (u in Rx) and (v in Ry)
        # 2. (u in x) and (v in Ry)
        # 3. (u in Rx) and (v in y)
        edges1 = _filter_edges(edges, xnodes, Ry)
        edges2 = _filter_edges(edges, Rx, ynodes)
        edges3 = _filter_edges(edges, Rx, Ry)
        edges1 = sorted(edges1, key=lambda x: x[2])
        edges2 = sorted(edges2, key=lambda x: x[2])
        edges3 = sorted(edges3, key=lambda x: x[2])
        w = min(edges1[0], edges2[0], edges3[0], key=lambda x: x[2])
        E.add(w)

        # Remove nodes from Rx, Ry
        Rx = Rx.difference({w[0], w[1]})
        Ry = Ry.difference({w[0], w[1]})

    # sum all of the weights
    weight = sum(list(zip(*E))[2])
    return E, weight

def _filter_edges(edges, x, y):
    """ Filters out edges that don't have nodes in x and y

    Parameters
    ----------
    edges : list, tuple
        List of edges (node1, node2, weight)
    x : set
        Set of nodes in sample x
    y : set
        Set of nodes in sample y

    Returns
    -------
    list, tuple
        List of tuples (node1, node2, weight)
    """
    new_edges = list()
    for e in edges:
        if ((e[0] in x and e[1] in y) or
            (e[1] in x and e[0] in y)):
            new_edges.append(e)
    return new_edges
