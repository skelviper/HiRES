# functions for HiRES project analysis
library(tidyverse)
library(Matrix)

# io functions
# read 4DN pairs file and return a dataframe
readPairs <- function(pairPath){
    pairs <- read_table(pairPath,col_names = FALSE, comment = "#") %>% select(1:5)
    names(pairs) <- c("readID","chrom1","pos1","chrom2","pos2")
    return(pairs)
}

# funcitons for CpG Clustering
readcolor2 <- function(color2_dir, chrList, splitname = ".pairs"){
    file_names <- dir(color2_dir, pattern = "*.color", recursive = F, full.names = T)
    # read color2 files and return a dataframe
    df <- read_table(file_names[1],col_names = F,col_types = cols(X1=col_character(),X2=col_double(),X3 = col_double())) %>% filter(X1 %in% chrList)
    df $ X1 <-as.character(df $ X1)
    names(df)<-c("chr","bin",strsplit(tail(strsplit(file_names[1],split="/")[[1]],1),split=splitname)[[1]][1])
    #print(df)
    # combine from second
    for (i in 2:length(file_names)) {
        if(i %% 100 == 0){
            print(paste0("Read ",i," cells done!"))
        }
        newdf <- read_table(file_names[i],,col_names = F,col_types = cols(X1=col_character(),X2=col_double(),X3 = col_double())) 
        if (dim(newdf)[1]<5){
            next
        }
        newdf <- newdf %>% filter(X1 %in% chrList)
        newdf $ X1 <-as.character(newdf $ X1)
        names(newdf)<-c("chr","bin",strsplit(tail(strsplit(file_names[i],split="/")[[1]],1),split=splitname)[[1]][1])
        df <- df %>% full_join(newdf,by=c("chr","bin"))
        }

    names(df)[1:2] <- c("chr","bin")
    return(df)
}

cleanCpGMatrix <- function(raw_cpg, bad_bin_threshold, bad_cell_threshold){
    NAs <- which(is.na(raw_cpg),arr.ind = T) %>% as_tibble() 
    rowToRemove <- NAs %>% group_by(row) %>% summarise(count = n()) %>% arrange(desc(count)) %>% filter(count>bad_bin_threshold) %>% pull(row)
    colToRemove <- NAs %>% group_by(col) %>% summarise(count = n()) %>% arrange(desc(count)) %>% filter(count>bad_cell_threshold) %>% pull(col)
    removeBad <- raw[-rowToRemove,-colToRemove]
    return(removeBad)
}
    
    
replaceCpGNA <- function(cleanCpGMatrix, cpgPath){
    cpg<- read_table2(cpgPath,col_names = c("chr","bin","cpg"))
    cpg$chr <- as.character(cpg$chr)
    
    toReplaceNA <- cleanCpGMatrix %>% left_join(cpg) %>% select(chr,bin,cpg,everything())
    replaceNAcoord <- which(is.na(toReplaceNA),arr.ind = T) %>% as_tibble() 
    replaceNA <- toReplaceNA
    for(i in 1:nrow(replaceNAcoord)){
        replaceNA[replaceNAcoord[i,1] %>% as.numeric,replaceNAcoord[i,2] %>% as.numeric] = replaceNA[replaceNAcoord[i,1] %>% as.numeric ,3]
    }
    replaceNA <- replaceNA %>% select(-cpg)
    return(replaceNA)
}

add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}