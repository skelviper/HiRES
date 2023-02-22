import time
import sys
import gzip
import os
import sys
import uuid
import json
import pandas as pd
from functools import partial
from pkgutil import get_data
from io import StringIO
from . import reference


'''
work with local pairs file generated from hickit
has 27 headline start with "#"
'''

# utils

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
## get sample name
def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:]) 

# parsers

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
def parse_gtf(filename:str) -> pd.DataFrame:
    # read gtf, get exons
    gencode = pd.read_table(filename, comment="#", header=None)
    gencode.columns="seqname source feature start end score strand frame group".split()
    return gencode
def parse_3dg(filename:str)->pd.DataFrame:
    # read in hickit 3dg file(or the .xyz file)
    # norm chr name alias

    ## get alias file in package
    ## reference is a "data module" with its own __init__.py
    dat = get_data(reference.__name__, "chrom_alias.csv")
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

# writers

def write_3dg(pairs:pd.DataFrame, outname:str):
    pairs.to_csv(outname, sep="\t",header=None)
    return 0
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

# [pipeline] json metrics recorder
def gen_record(record:dict, record_dir):
    uuid_str = uuid.uuid4().hex
    with open(os.path.join(record_dir, uuid_str+".json"), "w") as f:
        json.dump(record, f)