#! ~/miniconda3/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("usage: Rscript fraction.R input.tsv output.tsv", call.=FALSE)
}


library(tidyverse)
countMatrix<- read_table2(args[1])

if(dim(countMatrix)[2]<2){
    write_tsv(x=countMatrix,path=args[2])
} else {
    countMatrix <- countMatrix %>% mutate(overlapTimes = str_count(gene,",")+1) %>% select(gene,overlapTimes,everything())

    datNormed <-  lapply(countMatrix %>% select(-gene), function(x) x/countMatrix[[2]])
    datNormed <- cbind(countMatrix[1],datNormed) %>% select(-overlapTimes)
    datNormed <- datNormed %>% mutate(gene = strsplit(gene,",")) %>% unnest(gene) %>% group_by(gene) %>% summarise_each(funs(sum))

    write_tsv(x = datNormed,path = args[2])
}
