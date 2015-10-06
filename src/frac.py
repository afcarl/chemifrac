from __future__ import division
import numpy as np
import pandas as pd
import networkx as nx
import scipy


def frac(x, y, dm):
    """
    Parameters
    ----------
    x : np.array
       Vector of observations from first sample
    y : np.array
       Vector of observations from second sample
    dm : np.array
       Distance matrix

    Returns
    -------
    float
    A distance between graphs

    Notes
    -----
    The distance matrix needs to be in the same order as the vectors
    """
    # assert len(x) == len(y), 'vectors need to be the same length'
    # assert (x==0).sum()==0, 'x cannot have zeros'
    # assert (y==0).sum()==0, 'y cannot have zeros'

    D = len(x)
    d = 0
    for i in range(D-1):
        kv = ((x[i] - x[i+1:]) - (y[i] - y[i+1:]))**2
        dv = dm[i, i+1:]
        d += kv.dot(dv)
    return np.sqrt(d/D)


def sample_dm(mat, dm):
    """
    Computes the pairwise distance between all samples

    Parameters
    ----------
    mat : np.array
        OTU table
    dm : np.array
        OTU distance matrix

    Returns
    -------
    np.array:
        Sample distance matrix
    """
    r, c = mat.shape
    mat = np.log(mat)
    vamp_dm = np.zeros((r, r))
    for i in range(r):
        for j in range(i):
            vamp_dm[i, j] = frac(mat[i, :], mat[j, :], dm)
            print("i=%d, j=%d, dm[i,j]=%f"%(i, j, vamp_dm[i, j]))
        print(i)
    return vamp_dm


def connected_frac(x, y, dm, comp):
    """
    Parameters
    ----------
    x : np.array
       Vector of observations from first sample
    y : np.array
       Vector of observations from second sample
    dm : np.array
       Distance matrix
    comp : list
       List of connected components

    Returns
    -------
    float
    A distance between graphs

    Notes
    -----
    The distance matrix needs to be in the same order as the vectors
    """
    d = 0
    for u_ in range(len(comp)):
        u = comp[u_]
        for v_ in range(u_):
            v = comp[v_]
            kv = ((x[u] - x[v]) - (y[u] - y[v]))**2
            d += kv * dm[u, v]
    return d


def connected_dm(mat, dm):
    """
    Calculates connected components separately

    Parameters
    ----------
    mat : np.array
        OTU table
        rows = samples
        columns = features
    dm : np.array
       Distance matrix

    Returns
    -------
    float
    A distance between graphs

    Notes
    -----
    The distance matrix needs to be in the same order as the vectors
    """
    r, c = mat.shape
    mat = np.log(mat)
    D = c
    vamp_dm = np.zeros((r, r))
    # dmG = nx.Graph(dm)
    dmG = nx.from_scipy_sparse_matrix(scipy.sparse.dok_matrix(dm))
    connected_comps = list(nx.connected_components(dmG))
    for i in range(r):
        for j in range(i):
            d, D = 0, c
            for comp in connected_comps:
                # print(comp)
                vamp_dm[i, j] += connected_frac(mat[i ,:], mat[j, :],
                                                dm.values, list(comp))
                # print("i=%d, j=%d, dm[i,j]=%f"%(i, j, vamp_dm[i, j]))
            vamp_dm[i, j] = np.sqrt(vamp_dm[i, j]/D)
        # print(i)
    return vamp_dm
