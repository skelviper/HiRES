# calculate 3d procimity map from 3d genome structure
# .3dg file to .pairs file
# @author zliu

# load packages
from scipy import spatial
import numpy as np
import re
import argparse
import sys
import time
from numba import jit
import pandas as pd

parser = argparse.ArgumentParser(description='calculate 3d procimity map from 3d genome structure')
parser.add_argument('-i', '--input', help='input .3dg file', required=True,dest='tdgfilepath')
parser.add_argument('-o', '--output', help='output .pairs file', required=True,dest='pairsfilepath')
parser.add_argument('-d', '--distance', help='distance threshold', required=False, type=int, default=3,dest='distance')

args = parser.parse_args()

def read3dgfile(filepath:str):
    """
    read 3dg file and return pd.DataFrame
    """
    tdg = pd.read_csv(filepath,sep="\t",header = None)
    tdg.columns = ["chrom","pos","x","y","z"]
    tdg.sort_values(["chrom","pos"])

    return tdg
    
def sortPairs(pairs:pd.DataFrame):
    """
    sort upper triangle matrix 4DN pairs file
    """
    pairs.columns = ["readid","chrom1","pos1","chrom2","pos2","phase0","phase1"]
    pairs["pos1"] = pd.to_numeric( pairs["pos1"])
    pairs["pos2"] = pd.to_numeric( pairs["pos2"])
    pairs = pairs.sort_values(by=['chrom1', 'pos1', 'chrom2', 'pos2',"phase0","phase1"])
    pairs = pairs.drop_duplicates(subset=['chrom1', 'pos1', 'chrom2', 'pos2',"phase0","phase1"], keep='first')
    pairs = pairs.reset_index(drop=True)
    return pairs

start_time = time.time()

sys.stderr.write("Reading 3dg file\n")
tdg = read3dgfile(args.tdgfilepath)
nparray= np.array(tdg.values.tolist())

sys.stderr.write("Generating KD-tree\n")
kdtree = spatial.KDTree(nparray[:,2:])

sys.stderr.write("Querying pairs from KD-tree\n")
res = np.array(list(kdtree.query_pairs(r=args.distance)))

sys.stderr.write("Post processing of pairs file\n")

#@jit(nopython=True)
def arrayProcess(nparray,res):
    pos = nparray[:,:2]
    array = [np.append(pos[loc[0]], pos[loc[1]]) for loc in res]
    return array

df = pd.DataFrame(arrayProcess(nparray,res))
df.columns = ["chrom1","pos1","chrom2","pos2"]
df["phase0"] = df["chrom1"].apply(lambda x: 0 if re.search("pat",x) else 1)
df["phase1"] = df["chrom2"].apply(lambda x: 0 if re.search("pat",x) else 1)
df["chrom1"] = df["chrom1"].apply(lambda x: x[:-5])
df["chrom2"] = df["chrom2"].apply(lambda x: x[:-5])
df["readid"] = "."
df = df[["readid","chrom1","pos1","chrom2","pos2","phase0","phase1"]]

df = sortPairs(df)

sys.stderr.write("Writing outputs\n")
df.to_csv(args.pairsfilepath,sep="\t",index=False,header=None)

end_time = time.time()
sys.stderr.write("Duration: %ss\n" % (end_time - start_time))
