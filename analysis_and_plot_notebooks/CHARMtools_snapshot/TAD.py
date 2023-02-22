"""
Functions for single cell analysis in XingLab
@author zliu
@date: 2021/2/17
@version 0.0.2
"""

# Requirements,import some packages I'dont readlly understand why I need them
from __future__ import division, print_function
import scipy.stats as st
from scipy import ndimage
import cooler
import matplotlib.pyplot as plt
import lavaburst
import numpy as np
import pandas as pd
from hmmlearn import hmm

import seaborn as sns
sns.set_style('white')


def genome_coord2bed(genome_coord):
    """
    'chr1:35,900,000-40,900,000' to ['chr1',35900000,40900000]
    """
    return [genome_coord.split(sep=":")[0], int(genome_coord.replace(",", "").split(sep=":")[1].split(sep="-")[0]),
            int(genome_coord.replace(",", "").split(sep="-")[1])]


def getMatrixFromMCOOLs(filepath: "str", genome_coord, resolution=40000, balance=False):
    """
    intput: mcool file ,genome_coord(e.g. 'chr1:35,900,000-40,900,000'), resolution(should be included in mcoolfile)
    """
    cool = filepath+"::/resolutions/"+str(resolution)

    c = cooler.Cooler(cool)
    mat = c.matrix(balance=balance).fetch(genome_coord).astype("double")

    return mat


def lavaModularityProb(matrix, beta=10000, gamma=2, good_bin_nums=100):
    """
    description: calcualte modularity_score probability using lavaburst from single cell
    input: contact matrix
    output: list of modularity_score prob
    """
    good_bins = matrix.astype(bool).sum(axis=0) > good_bin_nums
    S = lavaburst.scoring.modularity_score(matrix, gamma, binmask=good_bins)
    model = lavaburst.SegModel(S)
    prob = model.boundary_marginals(beta)

    return prob


# convolution-randomwalk algorithm from schicluster
def gkern(kernlen=3, nsig=1):
    '''
    Returns a 2D Gaussian kernel array.
    '''

    interval = (2*nsig+1.)/(kernlen)
    x = np.linspace(-nsig-interval/2., nsig+interval/2., kernlen+1)
    kern1d = np.diff(st.norm.cdf(x))
    kernel_raw = np.sqrt(np.outer(kern1d, kern1d))
    kernel = kernel_raw/kernel_raw.sum()
    return kernel


def neighbor_ave_cpu(A, pad, mask="Mean", nsig=1):
    '''
    Convolution smooth
    '''

    if pad == 0:
        return A

    maskLength = pad*2 + 1
    maskKern = []
    if mask == "Mean":
        maskKern = np.array([1]*maskLength*maskLength).reshape(3, 3)
    if mask == "Gaussian":
        maskKern = gkern(maskLength, nsig)
    if mask == "Bilateral":
        pass

    return (ndimage.convolve(A, maskKern, mode='constant', cval=0.0) / float(maskLength * maskLength))


def random_walk_cpu(A, rp):
    """
    description:
    input:matrix,restart probability
    output:matrix after random walk
    """

    ngene, _ = A.shape
    A = A - np.diag(np.diag(A))
    A = A + np.diag(np.sum(A, axis=0) == 0)
    P = np.divide(A, np.sum(A, axis=0))
    Q = np.eye(ngene)
    I = np.eye(ngene)
    for i in range(30):
        Q_new = (1 - rp) * I + rp * np.dot(Q, P)
        delta = np.linalg.norm(Q - Q_new)
        Q = Q_new.copy()
        if delta < 1e-6:
            break
    return Q


def matImpute(A, pad, rp, mask):
    """
    description: run comvolution and random walk for a given matrix
    input: matrix_to_impute,convolution_pad_size,restart_probability_for_random_walk,convolution_mask_type
    """
    A = neighbor_ave_cpu(A, pad, mask)
    if rp == -1:
        Q = A[:]
    else:
        Q = random_walk_cpu(A, rp)

    return Q


def nice_ticks(ax, start, end, binsize, axis=(0, 1), tick_params=None):
    """
    I don't really know what this function for. (of course plot) zliu 1.24
    copy it from lavaburst example.
    """
    from matplotlib.ticker import MaxNLocator, FuncFormatter
    tick_locator = MaxNLocator(5)
    tick_formatter = FuncFormatter(
        lambda x, pos: '{:,}'.format(int(start + x*binsize)))
    if axis == 0 or axis == (0, 1):
        ax.xaxis.set_major_locator(tick_locator)
        ax.xaxis.set_major_formatter(tick_formatter)
    else:
        ax.set_xticks([])
    if axis == 1 or axis == (0, 1):
        ax.yaxis.set_major_locator(tick_locator)
        ax.yaxis.set_major_formatter(tick_formatter)
    else:
        ax.set_yticks([])

    if tick_params is not None:
        ax.tick_params(**tick_params)


def getModularityProbFromMCOOLs(mcoolpath, resolution, genome_coord, pad, rp, mask="Gaussian"):
    """
    for stat purposes, high modularity score means high insulation or low insulation score, shitty defs.

    output: return a list of modularity probobility of given genome coord and resolution.

    notes : if no convolution needed,set pad to 0, if no randomwalk needed,set rp to -1.
    """
    matrix = getMatrixFromMCOOLs(mcoolpath, genome_coord, resolution)
    matrix = matImpute(matrix, pad, rp, mask)
    prob = lavaModularityProb(matrix)
    return prob


def insulationScore(m, windowsize=500000, res=40000):
    """
    input: contact matrix,windowsize for sliding window, resolution of your contact matrix.
    ourput:

    """
    windowsize_bin = int(windowsize / res)
    score = np.ones((len(m)))
    for i in range(windowsize_bin, len(m) - windowsize_bin):
        with np.errstate(divide='ignore', invalid='ignore'):
            v = np.sum(m[max(0, i - windowsize_bin): i, i + 1: min(len(m) - 1, i + windowsize_bin + 1)]) / (np.sum(
                m[max(0, i - windowsize_bin):min(len(m) - 1, i + windowsize_bin + 1),
                    max(0, i - windowsize_bin):min(len(m) - 1, i + windowsize_bin + 1)]))
            if np.isnan(v):
                v = 1.0

        score[i] = v
    return score

def plotFromMCOOLs(mcoolpath, resolution, genome_coord, pad,rp,beta=10000,gamma=2,mask="Gaussian"):
    bedofcoord = genome_coord2bed(genome_coord)
    start = bedofcoord[1]
    end = bedofcoord[2]
    binsize = resolution

    matrix = getMatrixFromMCOOLs(mcoolpath, genome_coord, resolution)
    matrix = matImpute(matrix, pad, rp, mask)

    good_bins = matrix.astype(bool).sum(axis=0) > 10

    At = lavaburst.utils.tilt_heatmap(matrix, n_diags=500)

    S = lavaburst.scoring.modularity_score(matrix, gamma, binmask=good_bins)
    model = lavaburst.SegModel(S)

    prob = model.boundary_marginals(beta)

    f = plt.figure(figsize=(14, 10))
    ax = f.add_subplot(111)
    ax.matshow(np.log(At), cmap="Reds")
    ax.plot(np.arange(len(prob))-0.5, -10*prob, 'k', lw=1)
    ax.set_xlim([0, len(matrix)])
    ax.set_ylim([40, -10])
    ax.set_aspect(0.5)
    nice_ticks(ax, start, end, binsize, axis=0)

def QCbyContacts(celllist,minContacts=150000,resolution=40000):
    import cooler
    QCpassed = []
    for cell in celllist:
        c = cooler.Cooler(cell+"::/resolutions/"+str(resolution))
        if (c.info['sum']>minContacts):
            QCpassed.append(cell)
            
    return QCpassed

def calcInsulationCorr(listofmcool,genome_coord,resolution=40000,windowsize=200000):
    insulationScoreSingleCell = []
    for i in listofmcool:
        m = getMatrixFromMCOOLs(i,genome_coord,resolution)
        insulationScoreSingleCell.append(insulationScore(m,windowsize))
    # calculate pileup insulation score
    pileupMat = getMatrixFromMCOOLs(listofmcool[0],genome_coord,resolution)
    for i in listofmcool[1:]:
        pileupMat += getMatrixFromMCOOLs(i,genome_coord,resolution)

    insulationScoreBulk = insulationScore(pileupMat,windowsize)

    pearsonCorr = []
    for i in insulationScoreSingleCell:
        pearsonCorr.append(np.corrcoef(i, insulationScoreBulk)[0,1])
        
    return pearsonCorr,insulationScoreBulk,insulationScoreSingleCell

class pseudobulk:
    pearsonCorr = []
    insulationScoreBulk = []
    insulationScoreSingleCell = []
    hmmType = []
    resolution = 0
    cord = []
    mean = np.array([])
    std = np.array([])
 
    def __init__(self,listofmcool,genome_coord,resolution=40000,windowsize=200000):
        self.pearsonCorr,self.insulationScoreBulk,self.insulationScoreSingleCell = calcInsulationCorr(listofmcool,genome_coord,resolution=40000,windowsize=200000)
        self.resolution = resolution
        self.cord = [(i+1)*resolution for i in range(len(self.pearsonCorr))]
   
    def getPearson(self):
        return self.pearsonCorr
    def getInsulationBulk(self):
        return self.insulationScoreBulk
    def getInsulationSingleCell(self):
        return np.array(self.insulationScoreSingleCell)
    def getModel(self):
        return self.model
    def getHmmType(self,n_components = 4,covariance_type="full",n_iter = 1000):
        self.trainHMMmodel(n_components,covariance_type, n_iter)
        return self.hmmType

    def plotInsulationHeatmap(self):
        ax = sns.clustermap(self.insulationScoreSingleCell, vmin=0, vmax=0.3, col_cluster=False,figsize=(15,10))

    def trainHMMmodel(self,n_components,covariance_type, n_iter):
        hmmModel = hmm.GaussianHMM(n_components,covariance_type, n_iter)
        insulationMatrix = self.getInsulationSingleCell()
        self.std = np.std(insulationMatrix,axis = 0)
        self.mean = np.mean(insulationMatrix,axis = 0)
        self.model = hmmModel.fit(np.vstack((std,mean)).T)

        self.hmmType = self.model.predict(np.vstack((std,mean)).T)