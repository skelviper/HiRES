#! ~/miniconda3/bin/Rscript
# 
# dedup by R2-> extract readid -> filter pairend tabfastq
# @author zliu
#
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("usage: Rscript CTHiRES.extract_dedup_read.R 10 30 R2.bam.bed R1.tabfastq.gz R2.tabfastq.gz out.R1.tabfastq.gz out.R2.tabfastq.gz", call.=FALSE)
}

library(tidyverse)
options(bedtools.path="/share/home/zliu/miniconda3/envs/R4/bin/")
library(bedtoolsr)

#config
max_sep_distance = as.numeric(args[1])
min_mapping_qual = as.numeric(args[2])
R2bed_path = args[3]
R1.tabfastq_path = args[4]
R2.tabfastq_path = args[5]
R1.dedup.tabfastq_path = args[6]
R2.dedup.tabfastq_path = args[7]

R2bed <- read_table2(R2bed_path,col_names = F)

# when there is no such tag
if(dim(R2bed)[1]==0){
    R2bed %>% write_tsv(R1.dedup.tabfastq_path,col_names = F)
    R2bed %>% write_tsv(R2.dedup.tabfastq_path,col_names = F)
} else {
    R2bed <- R2bed %>% mutate(start = ifelse(X6=="+",X2,X3),end = start+1) %>% select(X1,start,end,X4,X5,X6) %>% arrange(X1,start,end)

    R2bed_cluster <- bedtoolsr::bt.cluster(R2bed,d=max_sep_distance)
    names(R2bed_cluster) <- c("chrom","start","end","readid","mapping_qual","strand","cluster")

    readid <- R2bed_cluster %>% filter(mapping_qual >= min_mapping_qual) %>% group_by(cluster) %>%
                    arrange(desc(mapping_qual)) %>%slice(1) %>% ungroup() %>% select(readid) %>% unique()

    R1.tabfastq <- read_table2(R1.tabfastq_path,col_names = c("readid","index","seq","qual")) %>% right_join(readid) %>% select(-index)
    R2.tabfastq <- read_table2(R2.tabfastq_path,col_names = c("readid","index","seq","qual")) %>% right_join(readid) %>% select(-index)

    R1.tabfastq %>% arrange(readid) %>% write_tsv(R1.dedup.tabfastq_path,col_names = F)
    R2.tabfastq %>% arrange(readid) %>% write_tsv(R2.dedup.tabfastq_path,col_names = F)
}

