#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Attempts to minimize \sum_{ij} A_{ij} |pi_i - pi_j| with eta-trick and spectral
ordering
"""
import warnings
import sys
import numpy as np
from scipy.sparse import issparse, coo_matrix, lil_matrix, find
from scipy.linalg import toeplitz

import matplotlib.pyplot as plt

from .spectral_embedding_ import spectral_embedding
from .utils import compute_score
from .spectral_ordering_ import SpectralBaseline


def plot_mat(X, title='', permut=None, true_pos=None):
    """
    Plot the permuted matrix X.
    Used to visualize the behavior of the iterates of
    spectral_eta_trick.
    """

    if permut is not None:
        if issparse(X):
            (iis, jjs, _) = find(X)
            invp = np.argsort(permut)
            pis = invp[iis]
            pjs = invp[jjs]
            # pis = permut[iis]
            # pjs = permut[jjs]

            # Xl = X.copy().tocsr()
            # Xl = Xl[permut, :]
            # Xl = Xl.T[permut, :].T
            # (pis, pjs, _) = find(Xl)
        else:
            Xl = X.copy()
            Xl = Xl[permut, :]
            Xl = Xl.T[permut, :].T

    fig = plt.figure(1)
    plt.gcf().clear()
    axes = fig.subplots(1, 2)
    if issparse(X):
        if permut is None:
            (pis, pjs, _) = find(X)
        axes[0].plot(pis, pjs, 'o', mfc='none')
    else:
        axes[0].matshow(Xl, interpolation='nearest')
    if permut is not None:
        if true_pos is None:
            axes[1].plot(np.arange(len(permut)), permut, 'o', mfc='none')
        else:
            true_perm = np.argsort(true_pos)
            true_inv_perm = np.argsort(true_perm)
            axes[1].plot(np.sort(true_inv_perm[permut]), true_pos[permut], 'o', mfc='none')
    plt.title(title)
    plt.draw()
    plt.pause(0.01)

    return


def spectral_eta_trick(X, n_iter=50, dh=1, score_function='1SUM',
                       dh_score=None, return_score=False,
                       do_plot=False, circular=False, norm_laplacian=None,
                       norm_adjacency=None, eigen_solver=None,
                       scale_embedding=False,
                       add_momentum=None,
                       true_pos=None):
    """
    Performs Spectral Eta-trick Algorithm from
    https://arxiv.org/pdf/1806.00664.pdf
    which calls several instances of the Spectral Ordering baseline (Atkins) to
    try to minimize 1-SUM or Huber-SUM (instead of 2-SUM)
    with the so-called eta-trick.
    """

    (n, n2) = X.shape
    assert(n == n2)

    if dh_score is None:
        dh_score = dh

    spectral_algo = SpectralBaseline(circular=circular,
                                     norm_laplacian=norm_laplacian,
                                     norm_adjacency=norm_adjacency,
                                     eigen_solver=eigen_solver,
                                     scale_embedding=scale_embedding)

    best_perm = np.random.permutation(n)
    best_score = compute_score(X, score_function=score_function, dh=dh_score, perm=best_perm, circular=circular)

    if issparse(X):
        if not isinstance(X, coo_matrix):
            X = coo_matrix(X)

        r, c, v = X.row, X.col, X.data
        eta_vec = np.ones(len(v))

        for it in range(n_iter):

            X_w = X.copy()
            X_w.data /= eta_vec

            if add_momentum:
                eta_old = eta_vec

            new_perm = spectral_algo.fit_transform(X_w)
            # new_perm = np.argsort(new_perm)

            if np.all(new_perm == best_perm):  # stopping criterion
                break

            p_inv = np.argsort(new_perm)
            if p_inv[0] > p_inv[-1]:  # convention to avoid alternating between one permutation and its flipped version
                p_inv = (n-1) - p_inv
                new_perm = np.argsort(p_inv)

            # if new_perm[0] > new_perm[-1]:  # convention to avoid alternating between one permutation and its flipped version
            #     new_perm = (n-1) - new_perm
            #     # new_perm *= -1
            #     # new_perm += (n-1)

            new_score = compute_score(X, score_function=score_function, dh=dh_score, perm=new_perm, circular=circular)
            if new_score < best_score:
                best_perm = new_perm  # keep best permutation so far

            eta_vec = abs(p_inv[r] - p_inv[c])
            # eta_vec **= 2  # !!!! THIS IS JUST A TEST FOR ROBUST 2 SUM !!!

            if circular:
                eta_vec = np.minimum(eta_vec, n - eta_vec)

            eta_vec = np.maximum(dh, eta_vec)  # solve Huber(dh) loss rather than 1-SUM

            if add_momentum:
                eta_vec = (1-add_momentum) * eta_vec + add_momentum * eta_old

            if do_plot:
                title = "it %d, score: %1.5e" % (it, new_score)
                plot_mat(X, permut=new_perm, title=title, true_pos=true_pos)

    else:
        eta_mat = np.ones((n, n))

        for it in range(n_iter):

            X_w = np.divide(X, eta_mat)

            new_perm = spectral_algo.fit_transform(X_w)

            # if new_perm[0] > new_perm[-1]:  # convention to avoid alternating between one permutation and its flipped version
            #     new_perm *= -1
            #     new_perm += (n-1)

            p_inv = np.argsort(new_perm)
            if p_inv[0] > p_inv[-1]:  # convention to avoid alternating between one permutation and its flipped version
                p_inv = (n-1) - p_inv
                new_perm = np.argsort(p_inv)

            if np.all(new_perm == best_perm):  # stopping criterion
                break

            new_score = compute_score(X, score_function=score_function, dh=dh_score, perm=new_perm, circular=circular)
            if new_score < best_score:
                best_perm = new_perm  # keep best permutation so far

            eta_mat = abs(np.tile(p_inv, n) - np.repeat(p_inv, n))

            if circular:
                eta_mat = np.minimum(eta_mat, n - eta_mat)

            eta_mat = np.reshape(eta_mat, (n, n))
            eta_mat = np.maximum(dh, eta_mat)  # solve Huber(dh) loss rather than 1-SUM

            if do_plot:
                title = "it %d, score: %1.5e" % (it, new_score)
                plot_mat(X, permut=new_perm, title=title)

    if return_score:
        return(best_perm, best_score)
    else:
        return(best_perm)


def spectral_eta_trick2(X, n_iter=50, dh=1, score_function='1SUM',                                 return_score=False,
                        do_plot=False, circular=False, norm_laplacian=None,
                        norm_adjacency=None, eigen_solver=None,
                        scale_embedding=False,
                        add_momentum=None,
                        avg_dim=3, avg_scaling=False,
                        true_pos=None,
                        dh_score=None):
    """

    THIS IS A MODIFIED EXPERIMENTAL VERSION OF THE ABOVE spectral_eta_trick FUNCTION.
    Instead of just using the Fiedler vector of A./eta to compute the next iterate eta_{t+1},
    it averages over several dimensions of the Laplacian embedding, in the spirit of
    https://arxiv.org/pdf/1807.07122.pdf.
    Preliminary experimental results indicate that this works better than the above code
    when circular=False only. Thus, we keep the official version from above in the SpectralEtaTrick
    method.

    Performs Spectral Eta-trick Algorithm from
    https://arxiv.org/pdf/1806.00664.pdf
    which calls several instances of the Spectral Ordering baseline (Atkins) to
    try to minimize 1-SUM or Huber-SUM (instead of 2-SUM)
    with the so-called eta-trick.

    Parameters
        ----------
        n_iter : int, default 50
            Number of iterations.

        score_function : string, default pSUM
            Which score we aim to minimize. Either '1SUM', '2SUM', 'Huber', 'R2S'
            (robust 2SUM function from the paper).
            If Huber or R2S, it is computer with the parameter dh provided.
            By design, the algorithm seeks to minimize the Huber loss. However,
            we keep the permutation that yields the best score amongst all, according
            to the score computed with score_function.
            
        dh : int, default 1
            Parameter for the Huber loss minimized.

        circular : boolean, default False
            Whether we wish to find a circular or a linear ordering.

        eigen_solver : string, default 'arpack'
            Solver for the eigenvectors computations. Can be 'arpack', 'amg', or
            'lopbcg'. 'amg' is faster for large sparse matrices but requires the
            pyamg package.

        add_momentum : Nonetype or float, default None.
            gamma parameter in Algorithm... from the paper.
            If gamma > 0, we set eta_{t+1} = gamma * eta_t + (1-gamma) * eta^*,
            where eta^* is the solution at iteration (t).

        avg_dim : int, default 1.
            Number of dimensions to use in the spectral embedding.
            If d = 1, it is the regular eta trick with eta = |pi_i - pi_j|.
            If d > 1, instead we sum |pi^k_i - pi^k_j| over the d first dimensions,
            where pi^k is the permutation that sorts the coordinates of the k-th dimension
            of the spectral embedding (not just the first, which is the Fiedler vector).
        
        avg_scaling : boolean, default True.
            If avg_dim > 1, the previous sum is weighted by the default scaling 1/(1+k)
            if avg_scaling = True.

        return_score : boolean, default False.
            Whether to return the best score (computed with score function) or not.
        
        norm_laplacian : string, default "unnormalized"
            type of normalization of the Laplacian. Can be "unnormalized",
            "random_walk", or "symmetric".

        norm_adjacency : str or bool, default 'coifman'
            If 'coifman', use the normalization of the similarity matrix,
            W = Dinv @ W @ Dinv, to account for non uniform sampling of points on
            a 1d manifold (from Lafon and Coifman's approximation of the Laplace
            Beltrami operator)
            Otherwise, leave the adjacency matrix as it is.
            TODO : also implement the 'sinkhorn' normalization

        scale_embedding : string or boolean, default True
            if scaled is False, the embedding is just the concatenation of the
            eigenvectors of the Laplacian, i.e., all dimensions have the same
            weight.
            if scaled is "CTD", the k-th dimension of the spectral embedding
            (k-th eigen-vector) is re-scaled by 1/sqrt(lambda_k), in relation
            with the commute-time-distance.
            If scaled is True or set to another string than "CTD", then the
            heuristic scaling 1/sqrt(k) is used instead.
        
    """

    (n, n2) = X.shape
    assert(n == n2)

    if n < 3:
        best_perm = np.arange(n)
        if return_score:
            return(best_perm, -1)
        else:
            return(best_perm)

    best_perm = np.arange(n)
    best_score = compute_score(X, score_function=score_function, dh=dh, perm=None)

    if issparse(X):
        if not isinstance(X, coo_matrix):
            X = coo_matrix(X)

        r, c, v = X.row, X.col, X.data
        eta_vec = np.ones(len(v))
        if add_momentum:
            eta_old = np.ones(len(v))

        for it in range(n_iter):

            X_w = X.copy()
            # eta_vec = np.ones(len(v))
            X_w.data /= eta_vec

            default_dim = 8
            if avg_dim > default_dim:
                default_dim = avg_dim + 1

            embedding = spectral_embedding(X_w, norm_laplacian=norm_laplacian,
                                           norm_adjacency=norm_adjacency,
                                           eigen_solver=eigen_solver,
                                           scale_embedding=scale_embedding,
                                           n_components=default_dim)

            new_perm = np.argsort(embedding[:, 0])

            # new_perm = spectral_algo.fit_transform(X_w)
            if np.all(new_perm == best_perm):
                break
            # if new_perm[0] > new_perm[-1]:
            #     embedding = embedding[::-1, :]
            #     new_perm *= -1
            #     new_perm += (n-1)

            new_score = compute_score(X, score_function=score_function, dh=dh, perm=new_perm)
            if new_score < best_score:
                best_perm = new_perm

            p_inv = np.argsort(new_perm)

            # eta_vec = abs(p_inv[r] - p_inv[c])
            eta_vec = np.zeros(len(r), dtype='float64')
            d_ = min(avg_dim, n-1)
            # eta_vec = np.sum(abs(embedding[r, :d_] - embedding[c, :d_]), axis=1)

            for dim in range(d_):
                # eta_mat = eta_mat + abs(np.tile(embedding[:, dim], n) - np.repeat(embedding[:, dim], n))
                d_perm = np.argsort(embedding[:, dim])
                d_perm = np.argsort(d_perm)
                eta_add = abs(d_perm[r] - d_perm[c])
                if circular:
                    eta_add = np.minimum(eta_add, n - eta_add)

                eta_add = np.maximum(dh, eta_add)

                if avg_scaling:
                    eta_add = eta_add * (1./(1 + dim))

                eta_vec += eta_add
                # eta_vec = eta_vec + eta_add
            #     eta_mat = eta_mat + abs(np.tile(d_perm, n) - np.repeat(d_perm, n))
            # eta_vec = np.sum(abs(embedding[r, :d_] - embedding[c, :d_]), axis=1)
            # if circular:
            #     # pass
            #     eta_vec = np.minimum(eta_vec, n - eta_vec)

            # pct = 50
            # this_delta = 50*n**(-3/2)
            # this_delta = 1.e-2
            # pct = 100 * ( 500. / n)
            # this_delta = np.percentile(eta_vec, pct)
            # this_delta = 1./10 * eta_vec.max()
            # eta_vec = np.maximum(this_delta, eta_vec)
            # print(eta_vec.std())
            # print(eta_vec.mean())
            # print(eta_vec.max())
            # eta_vec /= this_delta

            if do_plot:
                title = "it %d, score: %1.5e, delta=%1.3e" % (it, new_score, 0)
                plot_mat(X, permut=new_perm, title=title, true_pos=true_pos)

    else:
        eta_mat = np.ones((n, n))

        for it in range(n_iter):

            X_w = np.divide(X, eta_mat)

            default_dim = 8
            if avg_dim > default_dim:
                default_dim = avg_dim + 1

            embedding = spectral_embedding(X_w, norm_laplacian=norm_laplacian,
                                           norm_adjacency=norm_adjacency,
                                           eigen_solver=eigen_solver,
                                           scale_embedding=scale_embedding,
                                           n_components=default_dim)

            new_perm = np.argsort(embedding[:, 0])

            # new_perm = spectral_algo.fit_transform(X_w)
            # if new_perm[0] > new_perm[-1]:
            #     embedding = embedding[::-1, :]
            #     new_perm *= -1
            #     new_perm += (n-1)
            # if np.all(new_perm == best_perm):
            #     break

            new_score = compute_score(X, score_function=score_function, dh=dh, perm=new_perm)
            if new_score < best_score:
                best_perm = new_perm

            p_inv = np.argsort(new_perm)

            d_ = min(avg_dim, n-1)
            # eta_vec = np.sum(abs(embedding[r, :d_] - embedding[c, :d_]), axis=1)
            eta_mat = np.identity(n).flatten()
            for dim in range(d_):
                # eta_mat = eta_mat + abs(np.tile(embedding[:, dim], n) - np.repeat(embedding[:, dim], n))
                d_perm = np.argsort(embedding[:, dim])
                d_perm = np.argsort(d_perm)
                eta_add = abs(np.tile(d_perm, n) - np.repeat(d_perm, n))
                if circular:
                    eta_add = np.minimum(eta_add, n - eta_add)

                eta_add = np.maximum(dh, eta_add)

                if avg_scaling:
                    eta_add *= (1./(1 + dim))
                
                eta_mat = eta_mat + eta_add


            # eta_mat = abs(np.tile(p_inv, n) - np.repeat(p_inv, n))
            # if circular:
            #     # pass
            #     eta_mat = np.minimum(eta_mat, n - eta_mat)
            eta_mat = np.reshape(eta_mat, (n, n))
            # eta_mat = np.maximum(dh, eta_mat)

            if do_plot:
                title = "it %d, score: %1.5e" % (it, new_score)
                plot_mat(X, permut=new_perm, title=title)

    if return_score:
        return(best_perm, best_score)
    else:
        return(best_perm)
    

class SpectralEtaTrick():

    def __init__(self, n_iter=20, dh=1, return_score=False, circular=False,
                 norm_adjacency=None, eigen_solver=None, add_momentum=None,
                 do_plot=False, score_function='1SUM',
                 true_pos=None,
                 dh_score=None):
        self.n_iter = n_iter
        self.dh = dh
        self.return_score = return_score
        self.circular = circular
        self.norm_adjacency = norm_adjacency
        self.eigen_solver = eigen_solver
        self.add_momentum = add_momentum
        self.score_function = score_function
        self.dh_score = dh_score
        self.do_plot = do_plot
        self.true_pos = true_pos

    def fit(self, X):

        ordering_ = spectral_eta_trick(X, n_iter=self.n_iter, dh=self.dh,
                                       return_score=self.return_score,
                                       circular=self.circular,
                                       norm_adjacency=self.norm_adjacency,
                                       eigen_solver=self.eigen_solver,
                                       add_momentum=self.add_momentum,
                                       score_function=self.score_function,
                                       do_plot=self.do_plot,
                                       true_pos=self.true_pos,
                                       dh_score=self.dh_score)

        self.ordering = ordering_

        return self

    def fit_transform(self, X):

        self.fit(X)
        return self.ordering

