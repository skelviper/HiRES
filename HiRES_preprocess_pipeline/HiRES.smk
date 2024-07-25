####################################
#       HiRES pipeline             #
#@author Z Liu                     #
#@Ver 1.0                          #
#@date 2024/7/25                   #
####################################

#############CONFIG#################

import os

#input
SAMPLES = [i.split(sep='_')[0] for i in os.listdir("./Rawdata")]
#SAMPLES = os.listdir("./Rawdata")
#SAMPLES = ["d0-1"]

configfile: "HiRES_preprocess_pipeline/config.yaml"

#############RULE_ALL###############
"""
decide what you need for your down stream analysis.
"""
rule all:
    input:
        #preliminary split
        expand("processed/RNA_all/umibycell.{sample}.rna.R1.fq",sample=SAMPLES),
        #RNA part
        expand("result/RNA_Res/counts.{type}.{genome}.tsv",type=["gene","exon"],genome=["total","genome1","genome2"] if config["if_RNA_snp_split"] else ["total"]),
        expand("result/RNA_Res/counts.{type}.{genome}.format.tsv",type=["gene","exon"],genome=["total","genome1","genome2"] if config["if_RNA_snp_split"] else ["total"]),
        #Hi-C part pairs info
        expand("result/cleaned_pairs/c12/{sample}.pairs.gz",sample=SAMPLES),
        expand("result/dip_pairs/{sample}.dip.pairs.gz",sample=SAMPLES),
        #Hi-C part 3d info
        expand("processed/{sample}/3d_info/{sample}.{res}.align.rms.info",sample=SAMPLES if config["if_structure"] else [],res=["20k","50k","200k","1m"] if config["if_structure"] else []),
        expand("processed/{sample}/3d_info/{res}.{rep}.3dg", sample=SAMPLES if config["if_structure"] else [],
            res=["20k","50k","200k","1m"] if config["if_structure"] else [],
            rep=list(range(5)) if config["if_structure"] else []),
    threads: config["resources"]["generateStat_cpu_threads"] 
    shell:"""
        ./HiRES_preprocess_pipeline/HiRES_scripts/generateStat.sh
    """

    
############END_rule_all############


include: "rules/HiRES_split.rules"
include: "rules/scHiC_2dprocess.rules"
include: "rules/scHiC_3dprocess.rules"
include: "rules/HiRES_RNA.rules"


