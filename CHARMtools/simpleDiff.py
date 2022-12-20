# simpleDiff 
# Finding difference in two group of cooler

import re
from os import path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import diags

import gc

def getMatrixFromMCOOLs(filepath:str, genome_coord1:str,genome_coord2=None, resolution=40000, balance=False)->np.ndarray:
    """
    intput: mcool filepath ,
            genome_coord(e.g. 'chr1:35,900,000-40,900,000'), 
            resolution(should be included in mcoolfile)
    output: numpy 2d array
    """
    import cooler
    cool = filepath+"::/resolutions/"+str(resolution)

    c = cooler.Cooler(cool)
    if(genome_coord2 == None):
        genome_coord2 = genome_coord1
    matrix = c.matrix(balance=balance).fetch(genome_coord1,genome_coord2).astype("double")

    return matrix

def getMatrixFromCooler(filepath:str, genome_coord1:str,genome_coord2=None, resolution=40000, balance=False)->np.ndarray:
    """
    intput: cooler or mcool filepath, file type determined by last extension name.
            genome_coord(e.g. 'chr1:35,900,000-40,900,000'), 
            resolution(should be included in mcoolfile)
    output: numpy 2d array
    """
    import cooler
    
    if filepath.split(sep='.')[-1] == "mcool":
        return getMatrixFromMCOOLs(filepath, genome_coord1,genome_coord2, resolution, balance)
    
    c = cooler.Cooler(filepath)
    if(genome_coord2 == None):
        genome_coord2 = genome_coord1
    matrix = c.matrix(balance=balance).fetch(genome_coord1,genome_coord2).astype("double")

    return matrix

def getBandVec(mat:np.array,bins:int):
    return np.concatenate([np.diagonal(mat,offset = i) for i in range(bins)],axis=0)

def getBandMat(coolerPaths:list,genome_coord:str,bins=500):
    veclist = []
    for cellPath in coolerPaths:
        if path.exists(cellPath):
            mat = getMatrixFromCooler(cellPath,genome_coord1=genome_coord)
            veclist.append(getBandVec(mat,bins=bins))
    return np.array(veclist),mat.shape[0]

def multiple_testing_correction(pvalues, correction_type="FDR"):
    """
    Consistent with R - print
    correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05,
                                          0.069, 0.07, 0.071, 0.09, 0.1])
    from https://github.com/CoBiG2/cobig_misc_scripts/blob/master/FDR.py
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    sample_size = pvalues.shape[0]
    qvalues = empty(sample_size)
    if correction_type == "Bonferroni":
        # Bonferroni correction
        qvalues = sample_size * pvalues
    elif correction_type == "Bonferroni-Holm":
        # Bonferroni-Holm correction
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            qvalues[i] = (sample_size-rank) * pvalue
    elif correction_type == "FDR":
        # Benjamini-Hochberg, AKA - FDR test
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = sample_size - i
            pvalue, index = vals
            new_values.append((sample_size/rank) * pvalue)
        for i in range(0, int(sample_size)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            qvalues[index] = new_values[i]
    return qvalues

def simpleDiff(coolerPath1:list,coolerPath2:list,genome_coord:str,bins=500,resolution=20000,fdr_threshold=0.05):
    """
    simpleDiff main funciton, accept cooler path list, return bedpe format dataframe with statistics and FDR.
    """

    genome_coord_list = re.split("[-:]",genome_coord)
    mat1, reconstructMatrixShapeLength = getBandMat(coolerPath1[:],genome_coord = genome_coord,bins=bins)
    mat2, _ = getBandMat(coolerPath2[:],genome_coord = genome_coord,bins=bins)
    res = stats.ttest_ind(mat1,mat2)

    a= np.array(res)[1]
    statistics = np.array(res)[0]
    a[np.isnan(a)] = 1

    vec = multiple_testing_correction(a)
    diagonals = []
    diagonals_stats = []
    start=0
    for i in range(bins):
        diagonals.append(vec[start:start+reconstructMatrixShapeLength-i])
        diagonals_stats.append(statistics[start:start+reconstructMatrixShapeLength-i])
        start += reconstructMatrixShapeLength-i
    reMat = diags(np.array(diagonals),[i for i in range(bins)]).toarray()
    reMat_stats = diags(np.array(diagonals_stats),[i for i in range(bins)]).toarray()

    del mat1
    del mat2

    #format dataframe
    where = np.where((reMat < fdr_threshold) & (reMat > 0))
    where_list = (np.array(where)*resolution).tolist()
    where_list.append(reMat_stats[where].tolist())
    where_list.append(reMat[where].tolist())
    df = pd.DataFrame(where_list).T
    df.columns = ["start1","start2","stats","FDR"]
    df["end1"] = df["start1"] + resolution
    df["end2"] = df["start2"] + resolution
    df["chrom1"] = genome_coord_list[0]
    df["chrom2"] = genome_coord_list[0]
    df = df[["chrom1","start1","end1","chrom2","start2","end2","stats","FDR"]]
    
    del reMat
    gc.collect()
    
    return df

if __name__ == "__main__":
    import argparse
    import ray
    #import glob
    import time

    parser = argparse.ArgumentParser(description='Get simple diff between two cooler files')
    parser.add_argument('-i1', '--input1', nargs='+', help='cooler file path list for celltype1', required=True)
    parser.add_argument('-i2', '--input2', nargs='+', help='cooler file path list for celltype2', required=True)
    parser.add_argument('-b', '--bins', type=int, help='bins for matrix, e.g. for 20k resolution, 500bins needed for 10mb range', required=True)
    parser.add_argument('-t','--threads', help='number of threads', required=False, default=20)
    parser.add_argument('-r','--resolution', help='resolution of matrix, e.g. for 20k resolution, 20000', required=False, default=20000)
    parser.add_argument('-o', '--output', help='output file path', required=True)

    chroms = ["chr"+str(i) for i in range(1,19)]
    #chroms = ["chr19","chr18","chr17"]

    ray.init(num_cpus=int(parser.parse_args().threads))
    simpleDiff_remote = ray.remote(simpleDiff)

    start = time.time()
    futures = [simpleDiff_remote.remote(parser.parse_args().input1,parser.parse_args().input2,genome_coord=chrom,
    bins = parser.parse_args().bins,) for chrom in chroms]
    res = pd.concat(ray.get(futures), axis=0)
    res.to_csv(parser.parse_args().output,sep="\t",index=False)
    end = time.time()
    ray.shutdown()
    print("duration = ",end-start," seconds")
