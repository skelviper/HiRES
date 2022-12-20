# CHARMtools cell cycle analysis functions
# adapt from hic_basic from ychi
# @author zliu
# @date 2022-01-13

import numpy as np
import pandas as pd
from concurrent import futures
import gzip

from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import euclidean_distances, rbf_kernel

from .mdso.spectral_ordering_ import SpectralBaseline, SpectralOrdering
from python_tsp.heuristics.simulated_annealing import solve_tsp_simulated_annealing
from python_tsp.distances import euclidean_distance_matrix

from bisect import bisect_left, bisect_right
from itertools import dropwhile
from functools import partial

from .CHARMio import parse_pairs

def window_count(distances:pd.DataFrame, win_num:int,type="Nagano",resolution=100000)->pd.Series:
    # count distribution of distance array
    if type == "Nagano":
        ## breaks from Nagano2017:
        breaks = [0] + [1000*2**(0.125*i) for i in range(0, win_num)] 
    elif type == "Ps":
        ## ps curve-like breaks
        breaks = [i*resolution for i in range(win_num+1)]
    else:
        raise ValueError("type must be Nagano or Ps")

    window_count = []
    for index, break_ in enumerate(breaks):
        if index == win_num:
            count = len(distances[distances >= break_])
        else:
            count = len(distances[(distances >= break_) & (distances < breaks[index + 1])])
        window_count.append(count)
    window_count = pd.Series(window_count, 
        index = breaks)
    # normalized by all intra contacts
    return window_count/len(distances)

def dis_counts(cell_name:str,win_num=150,type="Nagano",resolution=100000):
    # work for 11 column table only
    # get cell's intra contact's distribution in Peter's window
    # using customized .pairs parser
    contacts = parse_pairs(cell_name)

    # get contact distance array
    intra = contacts.loc[contacts["chr1"] == contacts["chr2"]]
    distances = abs(intra["pos1"] - intra["pos2"])
    # count according to Peter's window range
    counts = window_count(distances,win_num=win_num,type=type,resolution=resolution)
    counts.name = cell_name
    #return counts
    return counts 

def euclid_kernel(cdps,n_components):
    # affinity based on PCA-euclid-distance
    # Input:
    #    n_components: number of PCA components used to calculate distance
    # Output:
    #    similarity matrix
    scaler = preprocessing.StandardScaler()
    scaled = scaler.fit_transform(cdps)
    pca = PCA()
    pca_res = pca.fit_transform(scaled)
    dm = euclidean_distances(pca_res[:,:n_components])
    sm = np.exp(-dm * 1/n_components)
    return sm

def euclid(cdps,n_components):
    scaler = preprocessing.StandardScaler()
    scaled = scaler.fit_transform(cdps)
    pca = PCA()
    pca_res = pca.fit_transform(scaled)
    dm = euclidean_distances(pca_res[:,:n_components])
    return dm

def c_rbf_kernel(cdps,n_components):
    scaler = preprocessing.StandardScaler()
    scaled = scaler.fit_transform(cdps)
    return rbf_kernel(scaled)

def _mdso(cdps:pd.DataFrame, annote:pd.DataFrame, n_components:int=6):
    """
    Input:
        cdps: contact decay profile, index must be sample name
        annote: annotaions, index must be sample name
    Output:
        list of sample_name as order
    """
    # calculate similarity matrix
    sm =pd.DataFrame(
        euclid_kernel(cdps.values,n_components),
        index = cdps.index,
        columns = cdps.index
        )
    so = SpectralOrdering()
    so_res = so.fit_transform(sm.values)
    return list(annote.iloc[so_res].index)


# start read here!

# funcitons for ordering 
def calc_cellcycle_ordering(filesp,startbin = None,endbin=None,method = "spectral",threads=24,n_components=6):
    """
    functions for ordering cellcycle
    Input:
        filesp : a pandas.Dataframe with a column called "pairs"
        method : "spectral" or "ra"(short for random annealing)
    Output:
        DataFrame with additional col "order_index"
    """
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(dis_counts,filesp["pairs"])
    ares = list(res)
    cdps = pd.DataFrame(ares)
    cdps.columns = cdps.columns.astype("string")
    cdps.index = filesp.index

    if startbin != None:
        subset = [0]+[i for i in range(startbin,endbin)]
        cdps = cdps.iloc[:,subset]

    if method == "spectral":
        order = _mdso(cdps,filesp,n_components)
        order_index = pd.Series(list(range(len(order))),index=order,name="order_index")
    elif method == "ra":
        dm = euclidean_distance_matrix(cdps.values)
        order = solve_tsp_simulated_annealing(dm)
        order_index = pd.Series(index=filesp.index[order[0]],data=list(range(filesp.shape[0])),name="order_index")
    else: 
        raise ValueError("method must be 'spectral' or 'ra'")

    ordered_filesp = pd.concat([filesp,order_index],axis=1)

    return ordered_filesp

# functions for calculating repli score
def _make_repli_dict(repli_chipf)->dict:
    # Parse and make per-chromosome sorted repli_score reference
    # Input:
    #    repli_chipf: file path
    # Output:
    #    {"chr1":{
    #        "starts" : [...],
    #        "ends" : [...],
    #        "values" : [...]
    #          }
    #    ...
    #    }
    repli_ref = pd.read_table(
        repli_chipf,
        header=None,
        names=["chrom","start","end","rep_value"]
    )
    repli_dict = {}
    for chrom, dat in repli_ref.groupby("chrom"):
        dat = dat.sort_values("start")
        repli_dict[chrom] = {
            "starts":dat["start"].values,
            "ends":dat["end"].values,
            "values":dat["rep_value"].values
        }
    return repli_dict
def _expand_mean(site,expand,starts,ends,values):
    lm = bisect_left(starts,site-expand)
    rm = bisect_right(ends,site+expand)
    return values[lm:rm].mean() if values[lm:rm].size > 0 else np.nan
def _repli_score(pairsf,repli_dict,expand=10_000,c1=1,p1=2,c2=3,p2=4):
    # Calculate replication score of Hi-C pairs file
    # Input:
    #    pairsf: file path
    #    repli_dict: reference segment-timing annote
    #    expand: length of expansion from query leg position
    #    c1,p1,c2,p2, key col index of pairs file
    # Output:
    #    repli_score: raw_repli_score defined by Nagano2017 
    #    annote_ratio: ratio of legs with valid reference timing
    try:
        with gzip.open(pairsf,"rt") as f:
            pos_rs = []
            for line in dropwhile(lambda x:x.startswith("#"),f):
                eles = line.strip().split("\t")
                pos_rs.append(_expand_mean(
                    int(eles[p1]),
                    expand,
                    repli_dict[eles[c1]]["starts"],
                    repli_dict[eles[c1]]["ends"],
                    repli_dict[eles[c1]]["values"]
                ))
                pos_rs.append(_expand_mean(
                    int(eles[p2]),
                    expand,
                    repli_dict[eles[c2]]["starts"],
                    repli_dict[eles[c2]]["ends"],
                    repli_dict[eles[c2]]["values"]
                ))
        pos_rs = pd.Series(pos_rs)
        repli_score = pos_rs[pos_rs > 0].count() / pos_rs.shape[0]
        annote_ratio = 1 - pos_rs[pos_rs.isna()].shape[0] / pos_rs.shape[0]
    except:
        return {"repli_score":np.NaN,"annote_ratio":np.NaN}
    else:
        return {"repli_score":repli_score,"annote_ratio":annote_ratio}


def calc_repli_score(filesp:pd.DataFrame,repli_chipf:str,expand=10_000,c1=1,p1=2,c2=3,p2=4,threads=24):
    """
    Calculate replication score of Hi-C pairs files
    Input:
        filesp: pd.DataFrame with columns "pairs"
        repli_chipf: file path of repli_chipf
        expand: length of expansion from query leg position
        c1,p1,c2,p2, key col index of pairs file
    Output:
        DataFrame with additional col "repli_score" and "annote_ratio"
    """
    repli_dict = _make_repli_dict(repli_chipf)
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(partial(_repli_score,repli_dict=repli_dict,expand=expand,c1=c1,p1=p1,c2=c2,p2=p2),filesp["pairs"])
    return pd.DataFrame(res)

#functions for calculate nearp and mitotic scores

def _contact_describe(pairsf:str,c1=1,p1=2,c2=3,p2=4) -> pd.Series:
    # get cell's basic statistics, defined in Nagano2017
    contacts = pd.read_table(pairsf, header=None, comment="#",low_memory=False)
    new_columns = list(contacts.columns)
    new_columns[c1] = "chr1"
    new_columns[p1] = "pos1"
    new_columns[c2] = "chr2"
    new_columns[p2] = "pos2"
    contacts.columns = new_columns
    intra = contacts.query(' chr1 == chr2 ')
    distances = abs(intra["pos1"] - intra["pos2"])
    
    all_ = len(distances[23_000 < distances])
    short = len(distances[(23_000 < distances) & (distances < 2_000_000)])
    mitotic = len(distances[(2_000_000 < distances) & (distances < 12_000_000)])
    farAvg = distances[(4_500_000 < distances) & (distances < 225_000_000)]

    mitotic_r = mitotic/all_
    short_r = short/all_
    
    # assign to different stages on Peter's cirtera
    if mitotic_r >= 0.3 and short_r <= 0.5:
        group = "Post-M"
    elif short_r > 0.5 and short_r + 1.8*mitotic_r > 1.0:
        group = "Pre-M"
    elif short_r <= 0.63:
        group = "G1"
    elif 0.63 < short_r <= 0.785:
        group = "early/mid-S"
    elif short_r > 0.785:
        group = "mid-S/G2"
    else:
        group = "blank"
    
    return {"near_p":short_r, "mitotic_p":mitotic_r, "farAvg":farAvg.mean(),"NaganoCellcycle":group }

def calc_premeta_score():
    pass

def calc_contact_describe(filesp:pd.DataFrame,c1=1,p1=2,c2=3,p2=4,threads=24):
    """
    Calculate cell's basic statistics, near% and mitotic% ,defined in Nagano2017
    Input:
        filesp: pd.DataFrame with columns "pairs"
        c1,p1,c2,p2, key col index of pairs file
    Output:
        DataFrame with additional col "short%", "mitotic%", "farAvg" and "group"
    """
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(partial(_contact_describe,c1=c1,p1=p1,c2=c2,p2=p2),filesp["pairs"])
    return pd.DataFrame(res)


def _QC_contacts_cell(pairsf,minContacts=50000):
    if(len(parse_pairs(pairsf)) > minContacts):
        return pairsf

def QCbycontacts(filesp,minContacts=50000,threads=36):
    with futures.ProcessPoolExecutor(threads) as pool:
        passname = pool.map(partial(_QC_contacts_cell,minContacts=minContacts),filesp["pairs"])
    res = filesp.query('pairs in @passname')
    return pd.DataFrame(res)