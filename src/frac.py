from __future__ import division
import numpy as np
import pandas as pd


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
            # print("i=%d, j=%d, dm[i,j]=%f"%(i, j, vamp_dm[i, j]))
        print(i)
    return vamp_dm
