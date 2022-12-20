load_mat <- function(filepaths,threads = 40){
    registerDoParallel(threads)
    mat <- foreach(filepath = filepaths ,.combine = "cbind",.errorhandling = "remove") %dopar%{
        example <- fread(filepath,header = FALSE,tmpdir = "./",col.names = c("chrom","pos1","chrom2","pos2","distance"))
        example <- example %>% mutate(band = pos2-pos1) %>% group_by(band) %>% 
            #mutate(distance_adj = scale(distance))
            #mutate(aveDistance = ave(distance), distance_adj = distance - aveDistanc)
            mutate(aveDistance = ave(distance), distance_adj = distance)
        example <- bintable %>% left_join(example) %>% select(chrom,pos1,pos2,distance_adj)
        example[is.na(example)] <- 10
        example <- cbind(example %>% head(dim(example)[1] / 2), example %>% tail(dim(example)[1] / 2) %>% select(distance_adj)) 
        names(example)[5] <- "distance_adj2"
        gc()
        example %>% select(distance_adj,distance_adj2)
    }
    return(as.matrix(mat))
}
convert_mat_to_contacts <- function(mat){
    mat[mat <= 3] <- 1
    mat[mat > 3] <- 0
    mode(mat) <- "integer"
    return(mat)
}
project1D <- function(bedpe){
    #accept a bedpe file file and project it to 1d-genome ,return bed file
    left <- bedpe%>% select(1:3)
    right <- bedpe %>% select(4:6)
    names(left) <- c("chrom","start","end")
    names(right) <- c("chrom","start","end")
    all <- rbind(left,right) %>%bt.sort() %>%  unique()
    return(all)
}
BiRD <- function(mat1,mat2,bintable,resolution=20000,bins=50,chunk_size=4000,alpha=0.05,low_thres=0.1,high_thres=0.9,threads=60){
    mat_names <- bintable %>% filter(str_detect(chrom,"mat")) %>% mutate(chrom = str_extract(chrom,"chr[0-9XY]+"),rowname = paste0(chrom,"_",as.character(pos1),"_",as.character(pos2))) %>% pull(rowname)
    passnames <- ((rowSums(mat1) / dim(mat1)[2]) > 0.1 ) & ((rowSums(mat1) / dim(mat1)[2]) < 0.9 ) & ((rowSums(mat2) / dim(mat2)[2]) > 0.1 ) & ((rowSums(mat2) / dim(mat2)[2]) < 0.9 )
    mat1 <- mat1[passnames,]
    mat2 <- mat2[passnames,]
    gc()
    split_num = floor(dim(mat1)[1] / chunk_size)
    alpha = 1- (1-alpha) ** (1/split_num)
    set.seed(42)
    registerDoParallel(threads)
    result <- foreach(i = seq(split_num), .combine = "c",.errorhandling = "pass") %dopar% 
    {   
        tempres <- list(SigDetDCF(t(mat1[((i-1)*chunk_size + 1):(i * chunk_size),]), t(mat2[((i-1)*chunk_size + 1):(i * chunk_size),]), foldlen = 4096, trunc = 0, MB = 1000, alpha = alpha, ReMax = 10, COMPU = 'Point', Pvm = 'Asy'))
        tempres[[1]]$index <- i
        gc()
        return(tempres)
    }
    tempData <- result
    # postprocess formating
    BSDCFres <- c()
    stpoint <- c()
    pvpoint <- c()
    for (i in seq(length(tempData))){
       if (class(tempData[[i]]) == "if" ){
           next
       }
       if (class(tempData[[i]][[1]][1]) == "character" ){
           next
       }
       if (names(tempData[[i]][[1]][1]) == "BSDCF_res"){
           tempData[[i]] <- tempData[[i]][[1]]
       }
       if (class(tempData[[i]][[1]][1]) != "character"){
           BSDCFres[[i]] <-  tempData[[i]]$BSDCF_res %>% rowwise %>% mutate(pos = list(seq(startind,endind) + chunk_size*(tempData[[i]]$index-1))) %>% select(-startind,-endind) %>% unnest(pos)
           pvpoint[[i]] <- tempData[[i]]$PvPoint %>% as.data.frame()
           stpoint[[i]] <- tempData[[i]]$StPoint %>% as.data.frame()
       }

    }

    differences <- cbind(BSDCFres %>% bind_rows() ,pvpoint %>% bind_rows() ,stpoint %>% bind_rows) # %>% select(3:5)
    names(differences) <- c("no1","no2","pos","pv","stat")
    differences <- differences %>% select("pos","pv","stat")
    location <- mat_names[passnames] %>% as.data.frame() %>% mutate(row_number())
    names(location) <- c("name","pos")
    diff_value <- rowSums(mat1) / dim(mat1)[2] - rowSums(mat2) / dim(mat2)[2]
    differences <- location %>% right_join(differences) %>% separate(name,into = c("chrom1","start1","start2")) %>% 
        mutate(start1=as.numeric(start1) - 1/2 * resolution,start2=as.numeric(start2) - 1/2 * resolution,chrom2 = chrom1,end1 = start1+ resolution,end2 = start2+resolution)
    diff_value <- diff_value %>% as.data.frame() %>% mutate(row_number())
    names(diff_value) <- c("diff","pos")

    differences <- differences %>% left_join(diff_value) %>% select(chrom1,start1,end1,chrom2,start2,end2,pv,diff)
                    
    return(differences)
}
