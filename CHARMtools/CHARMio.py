#global dependence
import numpy as np
import pandas as pd
import gzip
import os

#upsteram hic
def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:])
def parse_pairs(filename:str)->"Cell":
    '''
    read from 4DN's standard .pairs format
    compatible with all hickit originated pairs-like format 
    '''
    # read comments
    with gzip.open(filename,"rt") as f:
        comments = []
        chromosomes = []
        lengths = []
        for line in f.readlines():
            if line[0] != "#":
                break
            if line.startswith("#chromosome") or line.startswith("#chromsize"):
                chrom, length = line.split(":")[1].strip().split()
                chromosomes.append(chrom)
                lengths.append(int(length))
            if line.startswith("#columns:"):
                columns = line.split(":")[1].strip().split()
            ## comment lines are stored in dataframe.attrs["comment"]
            comments.append(line)
    dtype_array = {"readID":"category",
            "chr1":pd.CategoricalDtype(categories=chromosomes),
            "pos1":"int",
            "chr2":pd.CategoricalDtype(categories=chromosomes),
            "pos2":"int",
            "strand1":pd.CategoricalDtype(categories=["+","-"]),
            "strand2":pd.CategoricalDtype(categories=["+","-"]),
            "phase0":pd.CategoricalDtype(categories=["1","0","."]),
            "phase1":pd.CategoricalDtype(categories=["1","0","."]),
            "phase_prob00":"float",
            "phase_prob01":"float",
            "phase_prob10":"float",
            "phase_prob11":"float"}
    dtypes = {key:value for key, value in dtype_array.items() if key in columns}
    #read table format data
    pairs = pd.read_table(
        filename, 
        header=None, 
        comment="#",
        dtype=dtypes,
        names=columns
        )
    pairs.attrs["comments"] = comments
    pairs.attrs["name"], _ = divide_name(filename) # infer real sample name
    pairs.attrs["chromosomes"] = chromosomes
    pairs.attrs["lengths"] = lengths
    #assign column names
    #sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    return pairs

#Hi-C
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

def cooltoolsGetObsExp(filepath:str, genome_coord1:str,genome_coord2=None, resolution=40000, balance=False)->np.ndarray:
    """
    input: cooler or mcool path
    output: observe / expected matrix as np.ndarry 
    """
    import cooltools
    pass
    # if filepath.split(sep='.')[-1] == "mcool":
