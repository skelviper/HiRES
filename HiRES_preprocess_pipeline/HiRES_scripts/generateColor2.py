# Usage: generateColor2.py coolerPath cpg.bed output.tsv

import cooler
import numpy as np
import pandas as pd
import sys
sys.path.append('/shareb/zliu/analysis/hires_gastrulation/')
from CHARMtools import plot
from CHARMtools import compartment
from CHARMtools import loop
from CHARMtools import CHARMio

clr = cooler.Cooler(sys.argv[1])
linearCpGData = pd.read_csv(sys.argv[2],sep="\t",names=["chrom","start","end","cpg"],header=None)
chromlist = clr.chroms()[:][["name"]].values.T.tolist()[0][:19]

res = []
for chrom in chromlist:
    matrix = clr.matrix(balance=False).fetch(chrom).astype("int")
    #matrix = np.diag(np.ones(matrix.shape[0])) + matrix
    matrix[matrix>0] = 1
    linearCpGData_bychrom = linearCpGData.query("chrom == @chrom")
    linear_cpg_vector = linearCpGData_bychrom["cpg"].to_numpy()
    location = linearCpGData_bychrom[["chrom","start","end"]]
    location["scAB"] = np.dot(matrix,linearCpGData.query('chrom == @chrom')["cpg"].values) / (np.sum(matrix,axis=0)+1)
    res.append(location)
    
pd.concat(res).to_csv(sys.argv[3],sep="\t",header=None,index=None)