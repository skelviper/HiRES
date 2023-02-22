#global dependence
import numpy as np
import pandas as pd
import gzip
import os
from pkgutil import get_data
from io import StringIO
from . import ref
#import ref
from functools import partial

#upsteram hic
def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:])

def parse_pair_archieve(filename:str)->"Cell":
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

def parse_pairs(filename:str)->"Cell":
    '''
    read from 4DN's standard .pairs format
    compatible with all hickit originated pairs-like format 
    '''
    #comment lines are stored in dataframe.attrs["comment"]
    name_array = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1 phase_prob00 phase_prob01 phase_prob10 phase_prob11".split()
    #read comment line
    with gzip.open(filename,"rt") as f:
        comments = []
        for line in f.readlines():
            if line[0] != "#":
                break
            comments.append(line)
    #read table format data
    pairs = pd.read_table(filename, header=None, comment="#")
    pairs.attrs["comments"] = comments
    pairs.attrs["name"], _ = divide_name(filename) # get real sample name
    #assign column names
    pairs.columns = name_array[0:pairs.shape[1]]
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

def write_pairs(pairs:pd.DataFrame, out_name:str):
    '''
    write dataframe to tab delimited zipped file
    reserve comment lines, no dataframe index
    headers store in last comment line
    need to sort with upper triangle label
    '''
    #sys.stderr.write("write to %s\n" % out_name)
    with gzip.open(out_name,"wt") as f:
        pairs.attrs["comments"].pop()
        pairs.attrs["comments"].append("#columns:" + "\t".join(pairs.columns) + "\n")
        f.write("".join(pairs.attrs["comments"]))
        pairs.to_csv(f, sep="\t", header=False, index=False, mode="a")

def parse_gtf(filename:str) -> pd.DataFrame:
    # read gtf, get exons
    gencode = pd.read_table(filename, comment="#", header=None)
    gencode.columns="seqname source feature start end score strand frame group".split()
    return gencode

## norm name
def converter_template(c_in:str,ref_dict:pd.DataFrame):
    # a reat_table converter function
    #print(ref_dict)
    return ref_dict[c_in]
def fill_func_ref(template_func:callable, ref_file:str, index_col:str)->callable:
    # read in ref_file for template_fucn, generate new func
    # hope will boost new func's speed
    
    ## read in ref_file, get ref_dict in memory
    ref_df = pd.read_csv(ref_file, index_col=index_col)
    ref_dict = ref_df.iloc[:,0].to_dict()
    working_func = partial(template_func, ref_dict=ref_dict)
    return working_func

def parse_3dg(filename:str)->pd.DataFrame:
    # read in hickit 3dg file(or the .xyz file)
    # norm chr name alias

    ## get alias file in package
    ## reference is a "data module" with its own __init__.py
    dat = get_data(ref.__name__, "chrom_alias.csv")
    dat_f = StringIO(dat.decode())
    norm_chr = fill_func_ref(
                    converter_template,
                    dat_f,
                    "alias")
    ## read comments
    with open(filename) as f:
        comments = []
        for line in f.readlines():
            if line[0] != "#":
                break
            comments.append(line)
    ## read real positions
    s = pd.read_table(filename, 
                      comment="#",header=None,
                     index_col=[0,1],
                     converters={0:norm_chr})
    s.columns = "x y z".split()
    s.index.names = ["chr","pos"]
    ## assume last comment is backbone_unit
    if len(comments) > 0:
        s.attrs["comments"] = comments
        backbone_unit = float(comments[-1].split(":")[1].strip())
        s.attrs["backbone_unit"] = backbone_unit    
    return s

def write_3dg(pairs:pd.DataFrame, outname:str):
    pairs.to_csv(outname, sep="\t",header=None)
    return 0
