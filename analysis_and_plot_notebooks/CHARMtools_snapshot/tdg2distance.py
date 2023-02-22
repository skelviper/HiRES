# calculate 3d distacne from 3d genome structure
# .3dg file to .distance file
# @author zliu

# load packages
import numpy as np
import argparse
import pandas as pd
import tqdm

def read3dgfile(filepath:str):
    """
    read 3dg file and return pd.DataFrame
    """
    tdg = pd.read_csv(filepath,sep="\t",header = None)
    tdg.columns = ["chrom","pos","x","y","z"]
    tdg.sort_values(["chrom","pos"])

    return tdg

# args parse
parser = argparse.ArgumentParser(description="calculate 3d distance from 3d genome structure")
parser.add_argument("-i","--input",help="input 3dg file",required=True)
parser.add_argument("-o","--output",help="output distance file",required=True)
parser.add_argument("-b","--bins",help="bin size",type=int,default=50)

tdg = read3dgfile(parser.parse_args().input)

def calc_distance(tdg:pd.DataFrame,bins = 50):
    nparray= np.array(tdg.values.tolist(),dtype=np.str_)
    res = []
    for entry_index in tqdm.tqdm(range(len(nparray))):
        for i in range(bins):
            # only calculate distance intra chromosome
            if entry_index+i >= len(nparray) or nparray[entry_index][0] != nparray[entry_index+i][0]:
                break
            temp = list(nparray[entry_index][:2]) + list(nparray[entry_index + i][:2])
            distance = np.linalg.norm(nparray[entry_index][2:].astype(np.float64)-nparray[entry_index+i][2:].astype(np.float64))
            temp.append(distance)
            res.append(temp)
    return res

res = calc_distance(tdg,parser.parse_args().bins)
pd.DataFrame(res).to_csv(parser.parse_args().output ,sep="\t",index=False,header=None)


