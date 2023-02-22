import pandas as pd

from clean_isolated import clean_isolated
from CHARMio import parse_pairs, write_pairs
def cli(args):
    file_name, out_name1, out_name2, num_thread, up_dense, up_distance = \
        args.filename[0], args.output_file1, args.output_file2, int(args.num_thread), int(args.dense), int(args.distance)
    pairs = parse_pairs(file_name)
    dip_frame, hickit_frame = sep_clean(pairs, num_thread, up_dense, up_distance)
    write_pairs(dip_frame, out_name1)
    write_pairs(hickit_frame, out_name2)

# --------- working module --------
hap_word = {"1":"(mat)", "0":"(pat)"} #1 for maternal, 0 for paternal
def add_ap(data:pd.DataFrame, row_picker:pd.Series, col_picker:list, ap:str):
    # add appendix for any data subset, defined by raw and col picker
    # do inplace
    # in pandas, assign value must be done in single step, don't do chained indexing
    ap = [ap for i in range(0, len(data.loc[row_picker, :]))] # try to figure out a more elegant way
    data.loc[row_picker, col_picker].astype("string")
    data.loc[row_picker, col_picker] = data.loc[row_picker, col_picker].str.cat(ap)
def rm_hap(data:pd.DataFrame):
    data["chr1"] = data["chr1"].astype("string").str.extract(r'(chr[\dXY]+)')
    data["chr2"] = data["chr2"].astype("string").str.extract(r'(chr[\dXY]+)')
    return data
def sep_clean(pairs: pd.DataFrame, num_thread:int, up_dense:int, up_distance:int) -> tuple:
    # sep pairs haplotype and do clean isolated again
    # generate 9-column .pairs file for 2D analysis as well as longer 13-column .pairs for hickit 3d building
    # call the former .dip.pairs for simplicity
    dip_frame = pd.DataFrame()
    hickit_frame = pd.DataFrame()
    # blank DataFrame can concat any DataFrame
    # dealing 4 phase_prob columns one by one and stack together
    for code in "00 01 10 11".split():
        # pick phased lines(one of phase_probs >= 0.75)
        row_picker = pairs["phase_prob" + code] >= 0.75
        # set phase0 and phas1 column according to phase_probs
        pairs.loc[row_picker, "phase0"] = code[0]
        pairs.loc[row_picker, "phase1"] = code[1]
        hickit_frame = pd.concat([hickit_frame, pairs.loc[row_picker, :]])
        # set chr1 and chr2 column according to phase_probs, hickit_frame don't need this
        add_ap(pairs, row_picker, "chr1", hap_word[code[0]])
        add_ap(pairs, row_picker, "chr2", hap_word[code[1]])
        dip_frame = pd.concat([dip_frame, pairs.loc[row_picker, :]])
    dip_frame.attrs = hickit_frame.attrs = pairs.attrs
    dip_frame = clean_isolated(dip_frame, num_thread, up_dense, up_distance) # clean isolated contacts for newly generated chromosomes
    hickit_frame = clean_isolated(hickit_frame, num_thread, up_dense, up_distance)
    
    return dip_frame, hickit_frame