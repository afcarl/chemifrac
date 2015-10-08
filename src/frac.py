from __future__ import division
import numpy as np
import pandas as pd
import networkx as nx
import scipy
from skbio.stats.composition import closure

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
        for j in range(i):
            d += ((x[i] - x[j]) - (y[i] - y[j]))**2 * dm[i,j]
            print('xi=%f xj=%f yi=%f yj=%f d=%f' % (x[i],x[j],y[i],y[j],dm[i,j]))
    return np.sqrt(d/D)


def sample_dm(mat, dm):
    """
    Computes the pairwise distance between all samples

    Parameters
    ----------
    mat : pd.DataFrame
        Abundance table
        rows = samples
        columns = features
    dm : array_like
        Distance matrix between features

    Returns
    -------
    np.array:
        Sample distance matrix
    """
    num_samps, num_feats = mat.shape
    mat[mat==0] = 1 # add pseudo count for zeros
    mat = np.log(mat.apply(closure, axis=1))
    vamp_dm = pd.DataFrame(np.zeros((num_samps, num_samps)),
                           index=mat.index,
                           columns=mat.index)
    for i in range(num_samps):
        for j in range(i):
            u, v = mat.index[i], mat.index[j]
            d = frac(mat.loc[u, :].values,
                     mat.loc[v, :].values,
                     dm)
            vamp_dm.loc[u, v] = d
            print('u=%s v=%s d=%f'%(str(u),str(v), d))
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
                                                dm, list(comp))
                # print("i=%d, j=%d, dm[i,j]=%f"%(i, j, vamp_dm[i, j]))
            vamp_dm[i, j] = np.sqrt(vamp_dm[i, j]/D)
        # print(i)
    return vamp_dm
