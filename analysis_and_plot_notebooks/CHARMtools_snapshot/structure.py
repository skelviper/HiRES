from . import CHARMio
import numpy as np
from os import path


def get3dProximityStackMatrix(listCoolPath:list,genome_coord1:str,genome_coord2=None,resolution=20000):
    """
    input: 
    output: percentage matrix
    """
    ori_len = len(listCoolPath)
    listCoolPath = [cool for cool in listCoolPath if path.exists(cool)]
    if len(listCoolPath) < ori_len:
        import warnings
        warnings.warn("Some cool files are not found, proceeding with only coolfiles exsit in the list.")
    ave = np.average([CHARMio.getMatrixFromCooler(cool,genome_coord1,genome_coord2,resolution) for cool in listCoolPath],axis=0) / 2
    return ave

#split numpy array into different length numpy arrays
def split_array(array,length):
    return [array[i:i+length] for i in range(0, len(array), length)]