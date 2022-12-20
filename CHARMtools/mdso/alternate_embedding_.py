#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
For additional experiments we want to try other embeddings than the
Laplacian embedding.
Given a similarity matrix on sparse or numpy format, it creates an
embedding with LLE or Isomap or MDS.
'''
import warnings
import numpy as np
from scipy.sparse import issparse, isspmatrix, coo_matrix, identity
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh, lobpcg, eigs
from .utils import check_similarity, compute_laplacian,\
                       _graph_is_connected, _set_diag
from sklearn.utils import check_random_state, check_array
from sklearn.utils.extmath import _deterministic_vector_sign_flip
from sklearn import manifold
# from sklearn.manifold import MDS, Isomap, LocallyLinearEmbedding, TSNE


# Classical MDS
def get_dist_mat(csgraph, return_diag=True):

    if csgraph.ndim != 2 or csgraph.shape[0] != csgraph.shape[1]:
        raise ValueError('csgraph must be a square matrix or array')

    n = csgraph.shape[0]
    
    if isspmatrix(csgraph):
        if csgraph.format in ('lil', 'dok'):
            m = csgraph.tocoo()
            needs_copy = False
        else:
            m = csgraph
            needs_copy = True
        
        if m.format == 'dia':
            m = m.copy()
        else:
            m = m.tocoo(copy=needs_copy)
        m.data *= -1
        m.data -= m.min()

        # Center the distance matrix and take opposite
        w = m.sum(axis=0).getA1() - m.diagonal()
        m.data *= -1
        m.data += (1/n) * w[m.row]
        m.data += (1/n) * w[m.col]
        m.data -= (1/n**2) * w.sum()
        m.data *= 1/2

    else:
        m = np.array(csgraph)

        m = np.max(m) - m
        # m *=-1
        # m -= m.min()

        # # Center the distance matrix
        J = np.identity(n) - (1./n) * np.ones((n, n))
        m = -0.5 * J @ m @ J

    if return_diag:
        w = np.zeros(n)
        return m, w
    else:
        return m


def classical_MDS_embedding(adjacency, n_components=8, eigen_solver='arpack',
                            random_state=None, eigen_tol=1e-15,
                            norm_laplacian=False, drop_first=True,
                            norm_adjacency=False, scale_embedding=False,
                            verb=0):

    adjacency = check_similarity(adjacency, normalize=norm_adjacency)

    try:
        from pyamg import smoothed_aggregation_solver
    except ImportError:
        if eigen_solver == "amg":
            warnings.warn("The eigen_solver was set to 'amg', but pyamg is "
                          "not available. Switching to 'arpack' instead")
            # raise ValueError("The eigen_solver was set to 'amg', but pyamg "
            #                  "is not available.")

    if eigen_solver is None:
        eigen_solver = 'arpack'
    elif eigen_solver not in ('arpack', 'lobpcg', 'amg'):
        raise ValueError("Unknown value for eigen_solver: '%s'."
                         "Should be 'amg', 'arpack', or 'lobpcg'"
                         % eigen_solver)

    random_state = check_random_state(random_state)

    n_nodes = adjacency.shape[0]
    # Whether to drop the first eigenvector
    if drop_first:
        n_components = n_components + 1
    if n_components > n_nodes:
        print(" n_components ({}) > ({}) n_nodes. setting \
              n_components=n_nodes".format(n_components, n_nodes))
        n_components = n_nodes

    if not _graph_is_connected(adjacency):
        warnings.warn("Graph is not fully connected, spectral embedding"
                      " may not work as expected.")

    if (eigen_solver == 'arpack' or eigen_solver != 'lobpcg' and
       (not isspmatrix(adjacency) or n_nodes < 5 * n_components)):
        try:

            dist_mat, dd = get_dist_mat(adjacency,
                                        return_diag=True)

            (lambdas, diffusion_map) = eigsh(dist_mat, n_components, which='LA',
                                             tol=eigen_tol)
            idx = np.array(-lambdas).argsort()
            d = lambdas[idx]
            embedding = diffusion_map.T[idx]
            if scale_embedding:
                if scale_embedding == 'CTD':
                    embedding[1:] = (embedding[1:, :].T * np.sqrt(1./d[1:])).T
                    # embedding = embedding.T
                elif scale_embedding == 'heuristic':
                    embedding = embedding.T * np.sqrt(1./np.arange(
                        1, n_components+1))
                    embedding = embedding.T
                else:
                    embedding *= dd

    
        except RuntimeError:
            warnings.warn("arpack did not converge. trying lobpcg instead."
                          " scale_embedding set to default.")
            # When submatrices are exactly singular, an LU decomposition
            # in arpack fails. We fallback to lobpcg
            eigen_solver = "lobpcg"

    else:
        raise ValueError("So far, only eigen_solver='arpack' is implemented.")

    embedding = _deterministic_vector_sign_flip(embedding)

    if drop_first:
        return embedding[1:n_components].T
    else:
        return embedding[:n_components].T


# metric MDS

def metric_MDS_embedding(adjacency, n_components=8, eigen_solver=None,
                            random_state=None, eigen_tol=1e-15,
                            norm_laplacian=False, drop_first=True,
                            norm_adjacency=False, scale_embedding=False,
                            verb=0,
                            metric=True):

    adjacency = check_similarity(adjacency, normalize=norm_adjacency)

    dist_mat = adjacency
    if isspmatrix(dist_mat):
        dist_mat = dist_mat.toarray()
    dist_mat *= -1
    dist_mat -= dist_mat.min()

    method = manifold.MDS(n_components=n_components, metric=metric,
                          random_state=random_state,
                          dissimilarity='precomputed',
                          eps=0.0001,
                          max_iter=500)

    method.fit(dist_mat)
    embedding = method.embedding_

    return embedding

# TSNE

def tSNE_embedding(adjacency, n_components=8, eigen_solver=None,
                   random_state=None, eigen_tol=1e-15,
                   norm_laplacian=False, drop_first=True,
                   norm_adjacency=False, scale_embedding=False,
                   verb=0):

    adjacency = check_similarity(adjacency, normalize=norm_adjacency)

    dist_mat = adjacency
    if isspmatrix(dist_mat):
        dist_mat = dist_mat.toarray()
    dist_mat *= -1
    dist_mat -= dist_mat.min()

    if n_components > 3:
        tsne_method = 'exact'
    else:
        tsne_method = 'barnes_hut'
    method = manifold.TSNE(n_components=n_components, metric='precomputed',
                           random_state=random_state,
                           method=tsne_method)

    method.fit(dist_mat)
    embedding = method.embedding_[:, :]

    return embedding

def get_embedding(adjacency, n_components=8, eigen_solver=None,
                   random_state=None, eigen_tol=1e-15,
                   norm_laplacian=False, drop_first=True,
                   norm_adjacency=False, scale_embedding=False,
                   verb=0,
                   method='spectral'):
    drop_first = False
    if method == 'cMDS':
        embedding = classical_MDS_embedding(adjacency,
                                            n_components=n_components,
                                            eigen_solver=eigen_solver,
                                            random_state=random_state,
                                            eigen_tol=eigen_tol,
                                            norm_adjacency=norm_adjacency,
                                            norm_laplacian=norm_laplacian,
                                            drop_first=drop_first,
                                            scale_embedding=scale_embedding,
                                            verb=verb)

    elif method == 'MDS':
        embedding = metric_MDS_embedding(adjacency,
                                         n_components=n_components,
                                         eigen_solver=eigen_solver,
                                         random_state=random_state,
                                         eigen_tol=eigen_tol,
                                         norm_adjacency=norm_adjacency,
                                         norm_laplacian=norm_laplacian,
                                         drop_first=drop_first,
                                         scale_embedding=scale_embedding,
                                         verb=verb,
                                         metric=True)

    elif method == 'NMDS':
        embedding = metric_MDS_embedding(adjacency,
                                         n_components=n_components,
                                         eigen_solver=eigen_solver,
                                         random_state=random_state,
                                         eigen_tol=eigen_tol,
                                         norm_adjacency=norm_adjacency,
                                         norm_laplacian=norm_laplacian,
                                         drop_first=drop_first,
                                         scale_embedding=scale_embedding,
                                         verb=verb,
                                         metric=False)

    elif method == 'TSNE':
        embedding = tSNE_embedding(adjacency,
                                   n_components=n_components,
                                   eigen_solver=eigen_solver,
                                   random_state=random_state,
                                   eigen_tol=eigen_tol,
                                   norm_adjacency=norm_adjacency,
                                   norm_laplacian=norm_laplacian,
                                   drop_first=drop_first,
                                   scale_embedding=scale_embedding,
                                   verb=verb)
    
    else:
        drop_first = True
        embedding = spectral_embedding(adjacency,
                                       n_components=n_components,
                                       eigen_solver=eigen_solver,
                                       random_state=random_state,
                                       eigen_tol=eigen_tol,
                                       norm_adjacency=norm_adjacency,
                                       norm_laplacian=norm_laplacian,
                                       drop_first=drop_first,
                                       scale_embedding=scale_embedding,
                                       verb=verb)

    return embedding

if __name__ == '__main__':

    import os
    from time import time
    from scipy.sparse import coo_matrix
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mdso import SimilarityMatrix
    from mdso.utils import get_conn_comps
    from mdso.spectral_embedding_ import spectral_embedding

# ############ Synthetic example ################
    # Set parameters for data generation
    t0 = time()
    n = 500
    type_noise = 'gaussian'
    ampl_noise = 0.5
    type_similarity = 'CircularStrongDecrease'
    apply_perm = False
    # Build data matrix
    np.random.seed(44)
    data_gen = SimilarityMatrix()
    data_gen.gen_matrix(n, type_matrix=type_similarity, apply_perm=apply_perm,
                        noise_ampl=ampl_noise, law=type_noise)
    mat = data_gen.sim_matrix

    # mat = coo_matrix(mat)
    t1 = time()
    print("Generated similarity matrix -- {}".format(t1-t0))
    # Check it runs fine with different options
    norm_lap_opts = ['unnormalized', 'symmetric', 'random_walk']
    scaling_opts = [False, 'CTD', 'heuristic']
    norm_adj_opts = ['coifman', None]
    norm_adj_opts = [None]
    for norm_lap in [False]:
        for scale in scaling_opts:
            for norm_adj in norm_adj_opts:
                t_b = time()
                # embedding = classical_MDS_embedding(mat, norm_laplacian=norm_lap,
                #                                     scale_embedding=scale,
                #                                     norm_adjacency=norm_adj,
                #                                     drop_first=True,
                #                                     n_components=10,
                #                                     eigen_solver='arpack')

                # embedding = metric_MDS_embedding(mat, norm_laplacian=norm_lap,
                #                                     scale_embedding=scale,
                #                                     norm_adjacency=norm_adj,
                #                                     drop_first=True,
                #                                     n_components=3,
                #                                     eigen_solver='arpack')
                
                # embedding = spectral_embedding(mat, norm_laplacian=norm_lap,
                #                                     scale_embedding=scale,
                #                                     norm_adjacency=norm_adj,
                #                                     drop_first=True)          
                # 
                # 
                fig = plt.figure()
                axes = fig.subplots(2, 2)
                for idx, method in enumerate(['spectral', 'cMDS', 'MDS', 'TSNE']):
                    t_b = time()
                    if method == 'spectral':
                        drop_first = True
                    else:
                        drop_first = False
                    embedding = get_embedding(mat, norm_laplacian=norm_lap,
                                                scale_embedding=scale,
                                                norm_adjacency=norm_adj,
                                                drop_first=drop_first,
                                                n_components=5,
                                                eigen_solver='arpack',
                                                method=method)

                    print("Computed embedding with norm_lap : {}, \
                        scale_embedding : {} norm_adj : {}\
                        with method {} \
                        , in {}s.".format(norm_lap,
                                                                scale,
                                                                norm_adj,
                                                                method,
                                                                time()-t_b))
                    i_ax = idx // 2
                    j_ax = idx - (2) * i_ax
                    axes[i_ax, j_ax].scatter(embedding[:, 0], embedding[:, 1],
                                             c=np.arange(n))
                    axes[i_ax, j_ax].set_title(method)
                # fig = plt.figure()
                # ax = Axes3D(fig)
                # ax.scatter(embedding[:, 0], embedding[:, 2], embedding[:, 1],
                #         c=np.arange(n))
                # # plt.title("3d embedding of synthetic linear banded matrix")
                # plt.show()
                plt.tight_layout
                plt.show()

