#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tools for handling dense, sparse, connected and disconnected
similarity matrices
"""
import numpy as np
from scipy.sparse import issparse, coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components


def is_symmetric(m):
    """Check if a sparse matrix is symmetric
    (from Saullo Giovani)

    Parameters
    ----------
    m : array or sparse matrix
        A square matrix.

    Returns
    -------
    check : bool
        The check result.

    """
    if m.shape[0] != m.shape[1]:
        raise ValueError('m must be a square matrix')

    if not isinstance(m, coo_matrix):
        m = coo_matrix(m)

    r, c, v = m.row, m.col, m.data
    tril_no_diag = r > c
    triu_no_diag = c > r

    if triu_no_diag.sum() != tril_no_diag.sum():
        return False

    rl = r[tril_no_diag]
    cl = c[tril_no_diag]
    vl = v[tril_no_diag]
    ru = r[triu_no_diag]
    cu = c[triu_no_diag]
    vu = v[triu_no_diag]

    sortl = np.lexsort((cl, rl))
    sortu = np.lexsort((ru, cu))
    vl = vl[sortl]
    vu = vu[sortu]

    check = np.allclose(vl, vu)

    return check


def check_similarity(arr, normalize=False):
    """
    check that the matrix is square and symmetric, and
    normalize with Coifman's normalization (TODO : Sinkhorn-Knopp norm.).

    Parameters
    ----------
    arr : array or sparse matrix
        A square matrix with non-negative entries

    normalize : string or boolean or nonetype
        whether to perform normalization, and which type (so far, only Lafon
        and Coifman normalization [https://doi.org/10.1016/j.acha.2006.04.006]
        is implemented [TODO : sinkhorn knopp])

    Returns
    -------
    mat : bool
        A square matrix with non-negative entries, normalized with Lafon
        and Coifman's method if normalized='coifman'.
    """
    # check for squareness
    (n, m) = arr.shape
    if n != m:
        raise ValueError('the similarity matrix is not square')
    # check for symmetry
    if issparse(arr):
        if arr.format in ('lil', 'dok'):
            mat = arr.tocoo()
            needs_copy = False
        else:
            mat = arr
            needs_copy = True
        if (mat.data < 0).any():
            raise ValueError('the similarity matrix has negative entries')
        if not is_symmetric(mat):
            raise ValueError('specified similarity matrix is not\
                symmetric')
        if normalize:
            if not normalize == 'coifman':
                print("Warning: normalize argument is present but not Coifman.\
                      So far only Coifman's norm. is implemented!")
            w = mat.sum(axis=0).getA1() - mat.diagonal()
            mat = mat.tocoo(copy=needs_copy)
            isolated_node_mask = (w == 0)
            w = np.where(isolated_node_mask, 1, w)
            mat.data /= w[mat.row]
            mat.data /= w[mat.col]
            return(csr_matrix(mat, dtype='float64'))
        else:
            return(csr_matrix(mat, copy=needs_copy, dtype='float64'))

    else:
        if np.any(arr < 0):
            raise ValueError('the similarity matrix has negative entries')
        if not np.allclose(arr, arr.T, atol=1e-6):
            raise ValueError('specified similarity matrix is not\
                symmetric.')
        if normalize:
            if not normalize == 'coifman':
                print("Warning: normalize argument is present but not Coifman.\
                      So far only Coifman's norm. is implemented!")
            mat = np.array(arr)
            w = mat.sum(axis=0) - mat.diagonal()
            isolated_node_mask = (w == 0)
            w = np.where(isolated_node_mask, 1, w)
            mat /= w
            mat /= w[:, np.newaxis]
            return(mat)
        else:
            return(np.array(arr))


def get_conn_comps(mat, min_cc_len=1, return_n_cc=True):
    """
    Returns a list of connected components of the matrix mat by decreasing size
    of the connected components, for all cc of size larger or equal than
    min_cc_len
    """
    n_c, lbls = connected_components(mat)
    srt_lbls = np.sort(lbls)
    dif_lbls = np.append(np.array([1]), srt_lbls[1:] - srt_lbls[:-1])
    dif_lbls = np.append(dif_lbls, np.array([1]))
    switch_lbls = np.where(dif_lbls)[0]
    diff_switch = switch_lbls[1:] - switch_lbls[:-1]
    ord_ccs = np.argsort(-diff_switch)
    len_ccs = diff_switch[ord_ccs]
    ccs_l = []
    for (i, cc_idx) in enumerate(ord_ccs):
        if len_ccs[i] < min_cc_len:
            break
        ccs_l.append(np.where(lbls == cc_idx)[0])
    if return_n_cc:
        return ccs_l, n_c
    else:
        return ccs_l
