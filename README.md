# Code for HiRES analysis

<img src="label.png" width="50%" height="50%">

This folder project contains two folders, and I will describe their respective uses below.


## HiRES_preprocess_pipeline

The HiRES_preprocess_pipeline contains a snakemake pipeline for pre-processing HiRES data.

### Usage
1. Place the Rawdata folder and the HiRES_preprocess_pipeline folder in the same directory.
```bash
tree -h -L 2
# .
# ├── [ 252]  HiRES_preprocess_pipeline
# │   ├── [3.0K]  config.yaml
# │   ├── [  35]  envs
# │   ├── [ 248]  HiRES_scripts
# │   ├── [3.1K]  HiRES.smk
# │   ├── [1.0K]  LICENSE
# │   ├── [ 152]  README.md
# │   ├── [ 224]  rules
# │   ├── [ 361]  runHiRES_preprocess.sh
# │   └── [ 20K]  stat.ipynb
# └── [  33]  Rawdata
#     └── [  40]  OrgfE951001 -> ../../hires_new_test/Rawdata/OrgfE951001
```
2. Installation of environment
```bash
mamba create -n hires -c conda-forge -c bioconda python=3.8 snakemake=5.20.1 
```
3. Prepare files
HiRES_preprocess_pipeline relies on Hickit in addition to CHARMtools and the scripts that can be installed automatically, and you need to build the reference genome you need.
```bash
cd HiRES_preprocess_pipeline
vim config.yaml
```


## CHARMtools
The installation environment, CHARMtools, is a python package with some R code stored in the Rlib folder.

## analysis_and_plot_notebooks
analysis_and_plot_notebooks 中存放了项目分析中所用的python/R notebooks。