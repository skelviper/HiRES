#!/bin/bash

#usage: ./runHiRES_preprocess.sh 

cd ../
mkdir -p slurm_log
snakemake --cluster 'sbatch -w node03  --qos=medium --output=slurm_log/slurm-%j.out --cpus-per-task={threads} -t 7-00:00 -J HiRES!' --jobs 100 --resources nodes=100 --rerun-incomplete -s ./HiRES_preprocess_pipeline/HiRES.smk --keep-going

mkdir -p ./analysis
cp HiRES_preprocess_pipeline/stat.ipynb ./analysis/stat.ipynb
