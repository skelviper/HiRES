#####################################
# 3dg aligner for HiRES             #
# adapted from dip-c by Longzhi Tan #
# @author zliu                      #
# @date 2021/06/07                  #
#####################################

import sys
import getopt
import copy
import rmsd
import numpy as np
import itertools

def align(argv):
    # default parameters
    output_prefix = None

    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "o:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: hires3dAligner [options] <in1.3dg> <in2.3dg> ...\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -o STR        output prefix [no output]\n")
        sys.stderr.write("Output:\n")
        sys.stderr.write("  tab-delimited: homolog, locus, RMSD\n")
        sys.stderr.write("  additionally with \"-o\": 3DG files aligned to each other\n")

        return 1
    for o, a in opts:
        if o == "-o":
            output_prefix = a            
            
    # load 3dg files
    input_data = []
    num_structures = len(args)
    if num_structures < 2:
        sys.stderr.write("[E::" + __name__ + "] at least 2 structures are required\n")
        return 1
    counter = 0
    for input_filename in args:
        sys.stderr.write("[M::" + __name__ + "] reading 3dg file " + str(counter) + ": " + input_filename + "\n")
        input_data.append({})
        for input_file_line in open(input_filename, "rb"):
            input_file_line_data = input_file_line.strip().split()
            input_data[-1][(input_file_line_data[0], int(input_file_line_data[1]))] = [float(input_file_line_data[2]),float(input_file_line_data[3]),float(input_file_line_data[4])]
        counter += 1
    
    # find common particles
    common_loci = set(input_data[0])
    for input_structure in input_data[1:]:
        common_loci = common_loci.intersection(set(input_structure))
    num_loci = len(common_loci)
    common_loci = list(common_loci)
    common_data = []
    for input_structure in input_data:
        common_data.append([])
        for common_locus in common_loci:
            common_data[-1].append(input_structure[common_locus])
    sys.stderr.write("[M::" + __name__ + "] found " + str(num_loci) + " common particles\n")
        
    # subtract centroid
    common_data = np.array(common_data)
    centroid_data = []
    for i in range(num_structures):
        common_data[i] = np.array(common_data[i])
        centroid_pos = rmsd.centroid(common_data[i])
        common_data[i] -= centroid_pos
        centroid_data.append(centroid_pos)
    sys.stderr.write("[M::" + __name__ + "] found centroids for " + str(num_structures) + " structures\n")
    

    # data structure for store deviation for each pair
    # list of deviation [[num1,num2,deviation]....]
    list_of_deviations = []
    # calculate pairwise deviation and rotate
    deviations = np.empty((num_loci, 0), float)
    for i in range(num_structures):
        for j in range(num_structures):
            if j == i:
                continue
            # mirror image if needed
            mirror_factor = 1.0
            if rmsd.kabsch_rmsd(common_data[i], common_data[j]) > rmsd.kabsch_rmsd(common_data[i], -1.0 * common_data[j]):
                mirror_factor = -1.0
            # calculate deviation
            rotation_matrix = rmsd.kabsch(mirror_factor * common_data[j], common_data[i])
            if j > i:
                deviation = np.linalg.norm(np.dot(mirror_factor * common_data[j], rotation_matrix) - common_data[i], axis = 1).T
                deviations = np.c_[deviations, deviation]
                sys.stderr.write("[M::" + __name__ + "] median deviation between file " + str(i) + " and file " + str(j) + ": " + str(np.median(deviation)) + "\n")
                
                list_of_deviations.append((str(i)+str(j),np.median(deviation)))
                #list_of_deviations.append([j,i,deviation])
            
            # rotate
            if output_prefix is not None:
                # rotate j to align with i
                sys.stderr.write("[M::" + __name__ + "] aligning file " + str(j) + " to file " + str(i) + "\n")
                aligned_filename = output_prefix + str(j) + "_to_" + str(i) + ".3dg"
                aligned_file = open(aligned_filename, "wb")
                for input_locus in input_data[j]:
                    aligned_pos = np.dot((np.array(input_data[j][input_locus]) - centroid_data[j]) * mirror_factor, rotation_matrix) + centroid_data[i]
                    aligned_file.write("\t".join([input_locus[0], str(input_locus[1]), str(aligned_pos[0]), str(aligned_pos[1]), str(aligned_pos[2])]) + "\n")
                aligned_file.close()

    # summarize rmsd and print
    rmsds = np.sqrt((deviations ** 2).mean(axis = 1))
    totalrmsd = np.sqrt((rmsds ** 2).mean(axis = 0))
    sys.stderr.write("[M::" + __name__ + "] RMS RMSD: " + str(totalrmsd) + "\n")
    sys.stderr.write("[M::" + __name__ + "] median RMSD: " + str(np.median(rmsds,axis = 0)) + "\n")
    sys.stderr.write("[M::" + __name__ + "] writing output\n")
    for i in range(num_loci):
        sys.stdout.write("\t".join(map(str, [common_loci[i][0], common_loci[i][1], rmsds[i]])) + "\n")

    # summarize top 3 structur11es: 
    dict_of_diviations = dict(list_of_deviations)
    #sys.stderr.write(str(dict_of_diviations))
    str_num_structures = [str(i) for i in range(num_structures)]
    combinations = [i for i in list(itertools.product(str_num_structures,str_num_structures,str_num_structures)) if i[0]<i[1] and i[1]< i[2]]

    ser_conbination = 0
    min_deviation = 99999
    min_conbination = -1
    for i in range(len(combinations)):
        current_min_deviation = dict_of_diviations["".join([combinations[i][0],combinations[i][1]])] + \
                                dict_of_diviations["".join([combinations[i][1],combinations[i][2]])] + \
                                dict_of_diviations["".join([combinations[i][0],combinations[i][2]])]
        if min_deviation >= current_min_deviation:
            min_deviation = current_min_deviation
            min_conbination = ser_conbination
        ser_conbination+=1

    sys.stderr.write("this conbination persent min deviation: \n")
    sys.stderr.write(str(combinations[min_conbination]))
    sys.stderr.write("\n\n")
    
    #conbination_deviation = np.array([deviations[int(combinations[min_conbination][0])],\
    #                            deviations[int(combinations[min_conbination][1])],\
    #                            deviations[int(combinations[min_conbination][2])]])
    x = int(combinations[min_conbination][0])
    y = int(combinations[min_conbination][1])
    z = int(combinations[min_conbination][2])
    col1 = int(x*num_structures-(x**2+x)/2+y-x-1)
    col2 = int(y*num_structures-(y**2+y)/2+z-y-1)
    col3 = int(x*num_structures-(x**2+x)/2+z-x-1)
    conbination_deviation = np.array([deviations[:,col1],\
                                    deviations[:,col2],\
                                    deviations[:,col3]])
    conbination_rmsds = np.sqrt((conbination_deviation **2).mean(axis=1))
    conbination_totalrmsd = np.sqrt((conbination_rmsds**2).mean(axis=0))
    sys.stderr.write("[M::" + __name__ + "] top3 RMS RMSD: " + str(conbination_totalrmsd) + "\n")
    sys.stderr.write("[M::" + __name__ + "] top3 median RMSD: " + str(np.median(conbination_rmsds,axis = 0)) + "\n")
    sys.stderr.write("Example Structure presented: "+args[int(combinations[min_conbination][0])] + "\n")
    sys.stderr.write("[M::" + __name__ + "] Alldone!!!\n")


    return 0


if __name__ == "__main__":
    align(sys.argv)