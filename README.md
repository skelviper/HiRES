# Code for HiRES data analysis

[![DOI](https://zenodo.org/badge/580488541.svg)](https://zenodo.org/badge/latestdoi/580488541)

<img src="label.png" width="50%" height="50%">

## Introduction 

**Hi**-C and **R**NA-seq **e**mployed **s**imultaneously (HiRES) allows us to efficiently detect both the transcriptome and the 3D genome organization in a single cell. 

## Related Publications

Linking genome structures to functions by simultaneous single-cell Hi-C and RNA-seq, Science 380, 1070–1076 (2023)

Zhiyuan Liu*; Yujie Chen*; Qimin Xia*; Menghan Liu; Heming Xu; Yi Chi; Yujing Deng; Dong Xing✉

+ Raw data: PRJNA907173
+ Processed data: GSE223917

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
mamba activate hires
mamba env update --file HiRES_preprocess_pipeline/envs/main_env.yaml
```
3. Prepare files

HiRES_preprocess_pipeline relies on hickit(https://github.com/lh3/hickit/) and CHARMtools(https://github.com/skelviper/CHARMtools) in addition to softwares that can be installed automatically. Also you need to build index for RNA/DNA seperately on your version of reference genome.

```bash
cd HiRES_preprocess_pipeline
vim config.yaml
```

4. Run the pipeline
```bash
cd HiRES_preprocess_pipeline; ./runHiRES_preprocess.sh
```

5. generate statistics

    see analysis/stat.ipynb

## analysis_and_plot_notebooks

analysis_and_plot_notebooks holds the python/R notebooks used in the analysis of the project.
