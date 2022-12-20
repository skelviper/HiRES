#loop pileup analysis
#requirements
import cooler
import cooltools
import pandas as pd
import numpy as np

#loop should first be called and save in a tsv-like file contain at least chrom-start-end for both loop anchors
def cooltoolsGetOEPileUp(clr:cooler.Cooler,flank:int,resolution:int,expected:pd.DataFrame,loopAnchor:pd.DataFrame,arms:pd.DataFrame,nthreads:int)->np.ndarray:
    """
    cooltools warpper for OE pileup
    """
    loopAnchor.loc[:, 'mid1'] = (loopAnchor['start1'] + loopAnchor['end1'])//2
    loopAnchor.loc[:, 'mid2'] = (loopAnchor['start2'] + loopAnchor['end2'])//2

    stack = cooltools.pileup(clr, loopAnchor, view_df=arms, expected_df=expected, flank=flank,nproc=nthreads)

    mtx = np.nanmean(stack, axis=2)


    return mtx,stack
