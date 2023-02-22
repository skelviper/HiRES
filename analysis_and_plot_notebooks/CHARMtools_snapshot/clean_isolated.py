from concurrent import futures
import sys
import time
import bisect
from functools import partial

import pandas as pd
import numpy as np

from CHARMio import parse_pairs, write_pairs

def L_half(contact1, contact2):
    '''
    caculate L 0.5 norm of contact pair
    in general pandas groupby give chr1 < chr2
    light version, assume two contacts has same chromosome order
    '''
    d_x, d_y = abs(contact1.pos1-contact2.pos1), abs(contact1.pos2-contact2.pos2)
    return (np.math.sqrt(d_x) + np.math.sqrt(d_y))**2
def is_isolate(contact, sorted_contacts:"dataframe", up_dense, up_distance)->bool:
    '''
    check if contact is isolated, work on one chromosome pair
    contacts sorted by pos1 without index
    if two contact pos1 Eu distance > up_distance then L-0.5 distance > up_distance
    light version
    '''
    proximity = 0
    index = bisect.bisect_right(sorted_contacts["pos1"], contact["pos1"])
    for con in sorted_contacts.iloc[index-1::-1].itertuples(index=False):
        if abs(con.pos1-contact.pos1) > up_distance or proximity >= up_dense+1:
            break
        if L_half(con,contact) <= up_distance:
            proximity += 1 
    for con in sorted_contacts.iloc[index:].itertuples(index=False):
        if abs(con.pos1-contact.pos1) > up_distance or proximity >= up_dense+1:
            break
        if L_half(con,contact) <= up_distance:
            proximity += 1 
    return proximity < up_dense + 1
def clean_contacts_in_pair(contacts:"dataframe", up_dense, up_distance)->"dataframe":
    '''
    #for multiplexing
    '''
    sorted_contacts = contacts.sort_values(by="pos1",axis=0,ignore_index=True)
    mask = contacts.apply(is_isolate, axis=1, sorted_contacts = sorted_contacts,
                          up_dense=up_dense, up_distance=up_distance)
    sys.stderr.write("(%s, %s): %d --> %d\n" %(contacts.iloc[0]["chr1"], contacts.iloc[0]["chr2"], len(contacts),len(contacts[~mask])) )
    return contacts[~mask]
def cli(args):
    filename, out_name, num_thread, up_dense, up_distance = \
        args.filename[0], args.output_file, args.thread, args.dense, args.distance
    #working_function = partial(clean_isolated, num_thread=num_thread, up_dense=up_dense, up_distance=up_distance)
    pairs = parse_pairs(filename)
    write_pairs(clean_isolated(pairs, num_thread, up_dense, up_distance), out_name)
def clean_isolated(pairs, num_thread, up_dense, up_distance):
    t0 = time.time()
    input_data = ( value for key, value in pairs.groupby(["chr1", "chr2"]))
    working_func = partial(clean_contacts_in_pair, up_dense=up_dense, up_distance=up_distance)
    with futures.ProcessPoolExecutor(num_thread) as executor:
        res = executor.map(working_func, input_data)
    cleaned = pd.concat(res, axis=0)
    cleaned.attrs = pairs.attrs # groupby don't keep attrs
    print("clean_isolated: %d contacts removed in %s" % (len(pairs)-len(cleaned), pairs.attrs["name"]))
    sys.stderr.write("clean_isolated: finished in %.2fs\n"%(time.time()-t0))
    return cleaned
