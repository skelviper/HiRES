from argparse import ArgumentError
import os
import pandas as pd
import numpy as np

from CHARMio import parse_pairs, parse_3dg, write_3dg

def get_legs(pairs:pd.DataFrame)->pd.DataFrame:
    # get all legs of contacts
    leg_list1 = pairs[["chr1", "pos1"]]
    leg_list2 = pairs[["chr2", "pos2"]]
    leg_list1.columns = leg_list2.columns = ["chr","pos"]
    legs = pd.concat([leg_list1, leg_list2], ignore_index=True)
    #return legs
    return legs.set_index("chr")
def point_near_count(ref_points:np.array, scatter_points:np.array, radius:float)->np.array:
    # count number of points within radius
    # return array same length with ref_points
    left = scatter_points - radius
    right = scatter_points + radius
    scatter_points_expand = \
        pd.IntervalIndex.from_arrays(left=left, right=right)
    
    hit_num = []
    for ref_point in ref_points:
        res = scatter_points_expand.contains(ref_point)
        hit_num.append(
            res.astype(int).sum()
        )
    return np.array(hit_num)
def particle_evidence(chrom_structure:pd.DataFrame, legs:pd.DataFrame, chrom:str, max_clean_distance:float)->pd.DataFrame:
    # count number of contacts support structure particle
    # must provide particles' chromosome name
    # can't be used in groupby.transform
    ## output:
    ##     DataFrame same length/index of input chrom_structure
    try:
        chrom_legs = legs.loc[chrom, "pos"]
    except KeyError:
        # whole chromosome loss
        point_counts = np.zeros_like(
            chrom_structure.index.values
            )
        return  pd.Series(
            index = chrom_structure.index,
            data = point_counts
            )
    point_counts = point_near_count(
        chrom_structure.index.get_level_values("pos"),
        legs.loc[chrom, "pos"],
        max_clean_distance
    )
    return pd.Series(
        index = chrom_structure.index,
        data = point_counts
        )

# command entry point
def cli(args):
    structure, pairs, out_filename, clean_quantile, max_clean_distance = \
        args.filename, args.ref_filename, args.output, args.quantile, args.distance
    good_structure = clean3(
        s_name = structure,
        con_name = pairs,
        clean_quantile = clean_quantile,
        max_clean_distance = max_clean_distance 
    )
    write_3dg(good_structure, out_filename)
    return 0
def clean3(s_name, con_name, clean_quantile, max_clean_distance):
    # clean poorly supported particles
    # compatible with hickit
    # rescale x,y,z positions first

    # load data
    
    ## rescale 3dg
    s = parse_3dg(s_name)
    backbone_unit = s.attrs["backbone_unit"]
    norm_factor = 1.0/backbone_unit
    s[["x","y","z"]] = s[["x","y","z"]] * norm_factor
    ## read pairs
    pairs = parse_pairs(con_name)
    
    # get legs from contacts
    legs = get_legs(pairs)
    # count
    count_frames = []
    for name, group in s.groupby("chr"):
        count_frames.append(
            particle_evidence(
                group,
                legs,
                name,
                max_clean_distance
            )
        )
    s_contact_count = pd.concat(count_frames, axis=0)
    # get quantile
    gate_contact_num = s_contact_count.quantile(clean_quantile)
    # get good
    good_index = s_contact_count[s_contact_count > gate_contact_num].index
    good_3dg = s.reindex(good_index)
    return good_3dg
    