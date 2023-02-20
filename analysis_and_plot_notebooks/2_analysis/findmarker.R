#! ~/miniconda3/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

library(Seurat)
library(tidyverse)
source("/shareb//zliu//analysis/hires_gastrulation/CHARMtools/Rlibs/cellcycle_ana.R")
source("/shareb//zliu//analysis/hires_gastrulation/CHARMtools/Rlibs/plotFun.R")
source("/shareb//zliu//analysis/hires_gastrulation/CHARMtools/Rlibs/d3D.R")

hires <- readRDS("../../hires_mouse_emb_dev.new.rds")
hires <- subset(hires,subset = rmsd_20k < 1.5 & celltype != "mitosis" & celltype != "ExE ectoderm" & celltype != "ExE endoderm" & cellcycle_threshold != "M" & celltype != "blood")

options(scipen = 999999)
bintable <- read_tsv("../../figure3_related/d3d_res_analysis/bintable.20k.2m.tsv") %>% filter(pos2-pos1 < 2000000)
bintable <- bintable %>% filter(str_detect(chrom,pattern = "chr[0-9]+"))
binnames <- bintable %>% head(dim(bintable)[1] / 2) %>% mutate(chrom = str_extract(chrom,pattern = "chr[0-9]+"),pos = paste0(chrom,"_",pos1,"_",pos2)) %>% pull(pos)

celltypes <- hires[[]] %>% filter(rmsd_20k < 1.5,celltype != "mitosis") %>% 
    group_by(celltype) %>% summarise(count = n()) %>% arrange(desc(count)) %>% pull(celltype)
allcells <- hires[[]] %>% filter(rmsd_20k < 1.5,celltype != "mitosis") %>% pull(cellname) %>% sort()

chrom_now = args[1]

bintable_temp <- bintable %>% filter(chrom %in% paste0(chrom_now,c("(mat)","(pat)")))
filepaths <- paste0("./distance_chromosome/",allcells ,".",chrom_now,".tsv.gz")
binnames_temp <- bintable_temp %>% head(dim(bintable_temp)[1] / 2) %>% mutate(chrom = str_extract(chrom,pattern = "chr[0-9]+"),pos = paste0(chrom,"_",pos1,"_",pos2)) %>% pull(pos)

mat <- load_mat(filepaths,bintable_temp,type="all",threads = 80) %>% t() %>% as.data.frame() %>% suppressMessages()
rawmat <- mat[c(TRUE,FALSE),]
mat <- mat[c(FALSE,TRUE),]

rownames(mat) <-  paste0(sort(rep(allcells,2)),c("mat","pat"))
names(mat) <- binnames_temp
rownames(rawmat) <-  paste0(sort(rep(allcells,2)),c("mat","pat"))
names(rawmat) <- binnames_temp

test_wrapper <- function(cellnames_type1,cellnames_type2){
    mat1_mean_all <- rawmat[paste0(sort(rep(cellnames_type1,2)),c("mat","pat")),] %>% colMeans(na.rm = TRUE)
    mat2_mean_all <- rawmat[paste0(sort(rep(cellnames_type2,2)),c("mat","pat")),] %>% colMeans(na.rm = TRUE)

    mean_distance_all <-  cbind(mat1_mean_all,mat2_mean_all) %>% na.omit() %>% as.data.frame() %>% 
                            rownames_to_column("pos") %>% rowwise %>% 
                            filter(min(mat1_mean_all,mat2_mean_all) < 5) %>% ungroup()

    sig <- d3D( mat[paste0(sort(rep(cellnames_type1,2)),c("mat","pat")),pull(mean_distance_all,pos)],
                mat[paste0(sort(rep(cellnames_type2,2)),c("mat","pat")),pull(mean_distance_all,pos)],
       pull(mean_distance_all,pos),threads = 20,filter_type="pv",filter_thres =2,test.method = "wilcox")

    sig <- cbind(sig,mean_distance_all)
    return(sig)
}

sig_combine <- c()
for (celltype_i in celltypes){
    cellnames_this_type <- hires[[]] %>% filter(rmsd_20k < 1.5,celltype != "mitosis") %>% filter(celltype %in% celltype_i) %>% pull(cellname)
    celltype_i <- str_replace_all(celltype_i," ","_")
    restcells <- allcells[-match(cellnames_this_type,allcells)]

    sig <- test_wrapper(cellnames_this_type,restcells)
    sig %>% write_tsv(paste0("di_all/",celltype_i,".",chrom_now,".tsv"))
    
    sig_combine <- rbind(sig_combine,sig %>% mutate(celltype = celltype_i)%>% filter(pv < 0.05))
    gc()
    print(paste0(celltype_i," chrom ",chrom_now," is done!"))
}

for (seednow in c(41,42,43,44)){
    set.seed(seednow)
    cellnames_this_type <- hires[[]]%>% filter(rmsd_20k < 1.5,celltype != "mitosis") %>% sample_n(400) %>% pull(cellname)
    restcells <- allcells[-match(cellnames_this_type,allcells)]
    
    sig <- test_wrapper(cellnames_this_type,restcells)
    sig %>% write_tsv(paste0("di_all/Random",seednow,".",chrom_now,".tsv"))
    sig_combine <- rbind(sig_combine,sig %>% mutate(celltype = celltype_i)%>% filter(pv < 0.05))
    gc()
    print(paste0(seednow," chrom ",chrom_now," is done!"))
}

for (neuron_names_sample_num in c(50,100,150,200,250,300,400)){
    set.seed(42)
    cellnames_this_type <- hires[[]] %>% filter(celltype %in% c("early neurons")) %>% sample_n(neuron_names_sample_num) %>% pull(cellname)
    restcells <- allcells[-match(cellnames_this_type,allcells)]
    
    sig <- test_wrapper(cellnames_this_type,restcells)
    sig %>% write_tsv(paste0("di_all/Neuron_sample",neuron_names_sample_num,".",chrom_now,".tsv"))
    sig_combine <- rbind(sig_combine,sig %>% mutate(celltype = paste0("neuron",neuron_names_sample_num))%>% filter(pv < 0.05))
    gc()
    print(paste0("Neuronsample_",neuron_names_sample_num," chrom ",chrom_now," is done!"))
}

mat[,sig_combine%>% pull(pos) %>% unique()] %>% saveRDS(paste0("di_mat/",chrom_now,".mat.rds.gz"))
