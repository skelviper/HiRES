#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate various Similarity matrix
through the SimilarityMatrix methods
gen_matrix for synthetic data, and
gen_E_coli_matrix for DNA data.
"""
import numpy as np
from scipy.linalg import toeplitz


def gen_lambdas(type_matrix, n):
    '''
    Generates lambdas to define a toeplitz matrix with
    diagonal elements t_k = lambdas[k]
    There are four canonical matrices : LinearBanded, LinearStrongDecrease,
    CircularBanded, CircularStrongDecrease.
    '''
    array_lambdas = np.zeros(n)
    if type_matrix == 'LinearBanded':
        # Bandwidth = 10% ?
        cov = int(np.floor(n/10))
        array_lambdas[:cov] = cov - abs(np.arange(cov))

    elif type_matrix == 'LinearStrongDecrease':
        alpha = 0.1
        array_lambdas = np.exp(-alpha*np.arange(n))

    elif type_matrix == 'CircularBanded':
        # Bandwidth = 10% ?
        cov = int(np.floor(n/10))
        array_lambdas[:cov] = cov - abs(np.arange(cov))
        array_lambdas[-cov:] = array_lambdas[:cov][::-1]

    elif type_matrix == 'CircularStrongDecrease':
        alpha = 0.1
        array_lambdas = np.exp(-alpha*np.arange(n))
        p = int(np.floor(n/2))
        array_lambdas[-p:] = array_lambdas[:p][::-1]

    else:
        raise ValueError("Unrecognized type_matrix !")

    return(array_lambdas)


def gen_toeplitz_sim(lambdas):
    '''Build Toeplitz strong-R-matrix'''
    return(toeplitz(lambdas))


class SimilarityMatrix():

    # Apply permutation
    def apply_perm(self, perm):
        '''
        Apply a permutation to the similarity matrix.
        perm is given as a numpy array
        '''
        n_ = self.n
        # check size is ok
        if np.shape(perm)[0] != n_:
            raise ValueError('the size of the permutation matrix does not match that of the\
                             similarity matrix.')
        # check perm is a permutation
        if not (np.sort(perm) == np.arange(n_)).all():
            raise ValueError('perm is not considered as a'
                             'permutation matrix of [0; \cdots; n-1]')
        self.sim_matrix = self.sim_matrix[perm]
        self.sim_matrix = self.sim_matrix.T[perm]
        self.sim_matrix = self.sim_matrix.T
        return self

    # Add additive noise
    def add_sparse_noise(self, noise_prop, noise_eps,
                         law='uniform'):
        '''
        Adds a symetric sparse noise.
        noise_prop controls the support of the sparse noise
        noise_eps controls the eps amplitude of the noise
        '''

        n_ = self.n
        # first find a random support
        N = np.tril(np.random.rand(n_, n_))
        idx = np.where(N > noise_prop)
        N[idx] = 0
        # allocate value on the support
        (ii, jj) = np.where(N != 0)
        if law == 'gaussian':
            N[ii, jj] = noise_eps * np.abs(
                np.random.normal(0, 1, len(ii)))
        elif law == 'uniform':
            N[ii, jj] = noise_eps*np.random.rand(1, len(ii))
        # symetrize the noise
        N += N.T
        # Add noise to similarity matrix
        self.sim_matrix += N

        return self

    def gen_matrix(self, n, type_matrix='LinearBanded',
                   apply_perm=True, perm=None,
                   noise_prop=1, noise_ampl=0, law='uniform'):
        self.n = n
        lambdas = gen_lambdas(type_matrix, n)
        self.sim_matrix = gen_toeplitz_sim(lambdas)
        if apply_perm:
            if not perm:  # generate permutation if not provided by user
                perm = np.random.permutation(n)
            self.apply_perm(perm)
            self.true_perm = perm
        else:
            self.true_perm = np.arange(n)
        if noise_ampl > 0:
            normed_fro = np.sqrt(np.mean(self.sim_matrix**2))
            self.add_sparse_noise(noise_prop, noise_ampl*normed_fro, law=law)

        return self
