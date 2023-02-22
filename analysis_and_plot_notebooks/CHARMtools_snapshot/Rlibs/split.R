#! ~/miniconda3/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("usage: Rscript split.R cellname", call.=FALSE)
}


library(tidyverse)

cellname = args[1]
example <- read_tsv(paste0("distance2m/",cellname,".distance2m.gz"),col_names = FALSE)
names(example) <- c("chrom1","pos1","chrom2","pos2","distance")
example <- example %>% mutate(chrom = str_extract(chrom1,"chr[0-9]+"))
for (chrom_now in paste0("chr",seq(1:19))){
    print(chrom_now)
    example %>% filter(chrom == chrom_now) %>% select(- chrom) %>% write_tsv(paste0("distance2m_chromosome/",cellname,".",chrom_now,".tsv.gz"),col_names = FALSE)
}
