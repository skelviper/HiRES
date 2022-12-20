#!/bin/bash

#usage: ./runCHARM.sh

cd ../
mkdir -p slurm_log
snakemake --use-conda --cluster 'sbatch --exclude node03 --qos=medium --output=slurm_log/slurm-%j.out --cpus-per-task={threads} -t 7-00:00 -J CHARM!' --jobs 100 --resources nodes=100 --rerun-incomplete -s ./CHARM/CHARM.smk --keep-going

mkdir -p ./analysis
cp CHARM/stat.ipynb ./analysis/stat.ipynb
