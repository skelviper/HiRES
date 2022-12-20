#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
With an initial connected similarity matrix, our preprocessing might split
up the new similarity matrix is several components.
In the case of DNA reconstruction, the user might know that its data comes from
a single contigue. Hence it is of interest to merge the various contigue our
pre-processing might reveal.

Here twigs are supposed to be list of index forming a connected component of
the similarity_matrix.

The main function is merge
INPUT:
    - old_similarity
    - list connected components: list of index forming a connected component of
    the similarity_matrix.
OUTPUT:
    - a permutation of [1,...,n]
'''
import numpy as np
from scipy.sparse import issparse


def merge_conn_comp(l_connected, sim, h=15):
    if type(l_connected[-1]) == int:
        return(l_connected)
    n_cc = len(l_connected)
    l_out = l_connected.copy()
    a = np.zeros((2*n_cc, 2*n_cc))
    if issparse(sim):
        sim = sim.tolil()
    for i_cc in range(n_cc-1):
        for j_cc in range(i_cc+1, n_cc):
            for i_tip in ['begin', 'end']:
                for j_tip in ['begin', 'end']:
                    if i_tip == 'begin':
                        # filling beginning of i values
                        i_idxs = l_out[i_cc][:h]
                        this_i = 2*i_cc
                    elif i_tip == 'end':
                        # filling ending of i values
                        i_idxs = l_out[i_cc][-h:]
                        this_i = 2*i_cc + 1
                    if j_tip == 'begin':
                        # filling beginning of j values
                        j_idxs = l_out[j_cc][:h]
                        this_j = 2*j_cc
                    elif j_tip == 'end':
                        # filling ending of j values
                        j_idxs = l_out[j_cc][-h:]
                        this_j = 2*j_cc + 1

                    inter_sim = sim[i_idxs, :]
                    inter_sim = inter_sim.T[j_idxs, :].T
                    a[this_i, this_j] = inter_sim.sum()

    for _ in range(n_cc):
        # Find argmax
        (i_max, j_max) = np.unravel_index(a.argmax(), a.shape)
        if a[i_max, j_max] == 0:
            break
        # assert(i_max < j_max)
        if not (i_max < j_max):  # transpose i and j if i > j
            i_max_c = i_max
            i_max = j_max
            j_max = i_max_c
        i_tip = 'begin' if (i_max % 2 == 0) else 'end'
        i_cc = i_max // 2
        j_tip = 'begin' if (j_max % 2 == 0) else 'end'
        j_cc = j_max // 2
        # Merge sequences i and j
        s1 = l_out[i_cc]
        s2 = l_out[j_cc]
        if i_tip == 'begin':
            s1 = s1[::-1]
        if j_tip == 'end':
            s2 = s2[::-1]
        l_out[i_cc] = s1 + s2
        l_out[j_cc] = []
        # Update matrix after merging i and j (i becomes the merged sequence)
        # There are four possible cases (i is first arrow, j is second arrow)
        # --> --> ; --> <-- ; <-- --> ; <-- <--
        # Remove edges between i_cc and j_cc
        these_i = [2*i_cc, 2*i_cc + 1]
        these_j = [2*j_cc, 2*j_cc + 1]
        for this_i in these_i:
            for this_j in these_j:
                a[this_i, this_j] = 0
                a[this_j, this_i] = 0
        # a[these_i, :][:, these_j] = 0
        # a[these_j, :][:, these_i] = 0
        a_old = a.copy()
        if i_tip == 'begin':  # begin of new i becomes end of old i
            this_i = 2*i_cc + 1
            this_i_new = 2*i_cc
            a[this_i_new, :] = a_old[this_i, :]
            a[:, this_i_new] = a_old[:, this_i]
        if j_tip == 'begin':  # end of new i becomes end of old j
            this_j = 2*j_cc + 1
            this_i_new = 2*i_cc + 1
            a[this_i_new, :] = a_old[this_j, :]
            a[:, this_i_new] = a_old[:, this_j]
        elif j_tip == 'end':  # end of new i becomes beginning of old j
            this_j = 2*j_cc
            this_i_new = 2*i_cc + 1
            a[this_i_new, :] = a_old[this_j, :]
            a[:, this_i_new] = a_old[:, this_j]
        # erase component j
        for this_j in [2*j_cc, 2*j_cc+1]:
            a[:, this_j] = 0
            a[this_j, :] = 0
        # a = np.triu(a, k=1)

    l_out = [el for el in l_out if len(el) > 1]
    if len(l_out) == 1:
        l_out = l_out[0]
    return(l_out)
