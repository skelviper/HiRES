#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Evaluate the permutation given a ground truth
'''
import numpy as np
from scipy.stats import kendalltau
from scipy.sparse import issparse, coo_matrix
from scipy.linalg import toeplitz


def kendall_circular(true_perm, order_perm):
    '''
    TODO : make it faster for large n with a coarser grained slicing first,
    i.e., taking np.roll with a larger value than 1 and then zooming in.
    '''
    n = true_perm.shape[0]
    if (order_perm.shape[0] != n):
        print("wrong length of permutations in kendall_circular!")
    order_perm = true_perm[order_perm]
    id_perm = np.arange(n)
    scores = np.zeros(n)
    for i in range(n):
        scores[i] = abs(kendalltau(id_perm, order_perm)[0])
        order_perm = np.roll(order_perm, 1, axis=0)

    return(np.max(scores), np.argmax(scores))


def evaluate_ordering(perm, true_perm, criterion='kendall',
                      circular=False):
    '''
    evaluate the model.
    INPUT:
        - the ground truth permutation
        - the ordered_chain
    '''
    l1 = len(perm)
    l2 = len(true_perm)
    if not l1 == l2:
        print("Problem : perm of length {}, "
              "and true_perm of length {}".format(l1, l2))
        print("perm : {}".format(perm))
    if criterion == 'kendall':
        if circular:
            (score, _) = kendall_circular(true_perm, perm)
        else:
            score = abs(kendalltau(true_perm, np.argsort(perm))[0])
        return(score)


def compute_score(X, score_function='1SUM', dh=1, perm=None, circular=False):
    """ computes the p-sum score of X or X[perm, :][:, perm] if permutation
    provided
    """

    (n, _) = X.shape
    if issparse(X):
        if not isinstance(X, coo_matrix):
            X = coo_matrix(X)

        r, c, v = X.row, X.col, X.data

        if perm is not None:
            invp = np.argsort(perm)
            d2diag = abs(invp[r] - invp[c])
            # d2diag = abs(perm[r] - perm[c])
        else:
            d2diag = abs(r - c)

        if circular:
            d2diag = np.minimum(d2diag, n - d2diag)

        if not isinstance(dh, int):
            dh = int(dh)

        if score_function == '2SUM':
            d2diag **= 2
        elif score_function == 'Huber':
            is_in_band = (d2diag <= dh)
            # if circular:
            #     is_in_band += (d2diag >= n - dh)
            in_band = np.where(is_in_band)[0]
            out_band = np.where(~is_in_band)[0]
            d2diag[in_band] **= 2
            d2diag[out_band] *= 2 * dh
            d2diag[out_band] -= dh**2
        elif score_function == 'R2S':
            is_in_band = (d2diag <= dh)
            # if circular:
            #     is_in_band += (d2diag >= n - dh)
            in_band = np.where(is_in_band)[0]
            out_band = np.where(~is_in_band)[0]
            d2diag[in_band] **= 2
            d2diag[out_band] = dh**2

        prod = np.multiply(v, d2diag)
        score = np.sum(prod)

    else:
        if perm is not None:
            X_p = X.copy()[perm, :]
            X_p = X_p.T[perm, :].T
        else:
            X_p = X

        n = X_p.shape[0]
        d2diagv = np.arange(n)
        if circular:
            d2diagv = np.minimum(d2diagv, n - d2diagv)
        if score_function == '2SUM':
            d2diagv **= 2
        elif score_function == 'Huber':
            is_in_band = (d2diagv <= dh)
            # if circular:
            #     is_in_band += (d2diagv >= n - dh)
            in_band = np.where(is_in_band)[0]
            out_band = np.where(~is_in_band)[0]
            d2diagv[in_band] **= 2
            d2diagv[out_band] *= 2 * dh
            d2diagv[out_band] -= dh**2
        elif score_function == 'R2S':
            is_in_band = (d2diagv <= dh)
            # if circular:
            #     is_in_band += (d2diagv >= n - dh)
            in_band = np.where(is_in_band)[0]
            out_band = np.where(~is_in_band)[0]
            d2diagv[in_band] **= 2
            d2diagv[out_band] = dh**2

        D2diag_mat = toeplitz(d2diagv)
        prod = np.multiply(X_p, D2diag_mat)
        score = np.sum(prod)

    return score
