import sys
import time
import gzip
from concurrent import futures
from functools import partial
import pandas as pd


from CHARMio import parse_pairs, parse_gtf, write_pairs

def cli(args)->int:
    filename, gtf_file, out_name, thread = \
        args.filename[0], args.gtf_filename, args.out_name, args.num_thread
    # parsing ref gtf and pairs file
    pairs = parse_pairs(filename)
    # build in-mem exon index
    gtf = parse_gtf(gtf_file)
    ref = build_in_memory_index(get_exon(gtf))
    # do search
    cleaned = clean_splicing(pairs, ref, thread)
    write_pairs(cleaned, out_name)
def get_exon(gtf:pd.DataFrame) -> pd.DataFrame:
    # extract exon-gene_name from gtf table
    relevant = gtf.query('feature == "exon" & source == "HAVANA"') #using HAVANA only
    gene_id = relevant["group"].str.extract('gene_id "(ENSG[0-9]{11}.[0-9])";') #extract gene name from group
    gene_id.columns = ["gene_id"] # extract returns dataframe rather than series
    # don't mind strand
    return pd.concat([relevant.drop(["group","feature","source","score","strand","frame"],axis=1),gene_id],axis=1)
def build_in_memory_index(exons:pd.DataFrame) -> dict:
    # split by chr and using IntervalIndex to enable searching
    ref_dict = {key : value for key, value in exons.groupby("seqname")}
    # build index by chromosome
    for chromosome in ref_dict:
        # using start, end attrs as index
        by_chr_table = ref_dict[chromosome]
        bed_tuple = by_chr_table.set_index(['start','end']).index 
        bed_interval_index = pd.IntervalIndex.from_tuples(bed_tuple)
        by_chr_table.index = bed_interval_index
        ref_dict[chromosome] = by_chr_table.drop(["start","end","seqname"],axis=1)
    sys.stderr.write("hires_utils::clean_splicing: index done.\n")
    return ref_dict
def legs_co_gene(contact:pd.Series, chromosome:str, ref_dict:dict)->bool:
    # whether two legs of contacts are in same gene's exon
    # must be intra contacts
    pos1, pos2 = contact["pos1"], contact["pos2"]
    result = ref_dict[chromosome][ref_dict[chromosome].index.contains(pos1)]
    if len(result[result.index.contains(pos2)]) > 0:
        return True
    else:
        return False
def search_chromosome(contacts_at_chromosome:tuple, ref:dict)->pd.DataFrame:
    # search whole chromosome using F::legs_co_gene
    # single chr searching for multi-process calling
    chromosome, contacts = contacts_at_chromosome[0], contacts_at_chromosome[1]
    hit_index = contacts.apply(legs_co_gene, chromosome=chromosome, ref_dict=ref, axis=1)
    return contacts[hit_index]

def clean_splicing(pairs:pd.DataFrame, ref:dict, thread:int)->pd.DataFrame:
    '''
    clean contacts from splicing
    '''
    #intra = pairs.query('chr1 == chr2') # only search for intra
    intra = pairs.loc[pairs['chr1'] == pairs['chr2']]
    working_func = partial(search_chromosome, ref=ref) # pool.map can't take multiple iterable as arguments
    input_data = [(chromosome, per_chr_contacts) for chromosome, per_chr_contacts in intra.groupby("chr1")] # pool.map can't take additional keyword argument
    sys.stderr.write("hires_utils::clean_splicing: input parsed, search in %d thread\n" % thread)
    # do multi-thread searching
    with futures.ProcessPoolExecutor(thread) as pool:
        res = pool.map(working_func, input_data)
    result = pd.concat(res, axis=0) # contained contacts
    cleaned = pairs.drop(result.index, axis=0) # clean contacts
    print("clean_splicing: %d contacts removed in %s\n" %(len(result), pairs.attrs["name"]) )
    return cleaned

if __name__ == "__main__":
    # infile = "/shareb/zliu/project/202212/hires_pipe/result/cleaned_pairs/c12/OrgfE951001.pairs.gz"
    # gtf_file = "/share/Data/public/ref_genome/mouse_ref/M23/raw_data/annotation.gtf"
    # outfile = "/shareb/zliu/project/202212/hires_pipe/result/cleaned_pairs/c123/OrgfE951001.pairs.gz"
    infile = sys.argv[1]
    gtf_file = sys.argv[2]
    outfile = sys.argv[3]
    thread = 6

    pairs = parse_pairs(infile)
    # build in-mem exon index
    gtf = parse_gtf(gtf_file)
    ref = build_in_memory_index(get_exon(gtf))
    # do search
    cleaned = clean_splicing(pairs, ref, thread)
    write_pairs(cleaned, outfile)
