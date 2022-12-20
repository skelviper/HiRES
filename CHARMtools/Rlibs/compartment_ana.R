# Code for analysis of compartment dynamics
library(ggrepel)

# functions for discover scAB relationship with RNA expression (normally marker gene have higher scAB )

# Usage:
# genes <- read_table("~/share/Data/public/ref_genome/mouse_ref/M23/raw_data/gene.bed",col_names = c("chr","start","end","id","name","strand")) %>% group_by(name) %>% slice(1)
# resolution =1000000
# Idents(hires) <- hires$celltype
# options(repr.plot.width = 5,repr.plot.height=5,repr.plot.res = 200)
# meso_ecto_data <- getCpGDiffVocanoPlotData(hires,c("neural ectoderm"),c("early mesoderm"),genes,1000000,slot_cpg = "scab")
# meso_ecto_plot <- plotCpGDiffVocanoPlot(meso_ecto_data,"ectoderm vs. mesoderm",text_x_Pos = 20,scale_factor = 1e5)
# meso_ecto_plot+ theme_Publication()

getCpGDiffVocanoPlotData <- function(hires,ident1,ident2,genes,resolution=100000,slot_cpg = "cpg",slot_RNA = "SCT"){
    # calculate cpg diff
    DefaultAssay(hires) <- slot_cpg
    cpg_diff <- FindMarkers(hires,`ident.1` = ident1,`ident.2` = ident2,slot = "data",logfc.threshold = 0,mean.fxn = rowMeans,test.use = "wilcox")
    cpg_diff <- cpg_diff %>% rownames_to_column("region") %>% 
    mutate(chr = str_extract(region,"chr[0-9]+"),pos = as.numeric(str_extract(region,"[0-9.]+e.[0-9]{2}"))) 
    names(cpg_diff)[3] <- "mean_diff"
    # calculate rna diff
    DefaultAssay(hires) <- slot_RNA
    diff_SCT <- FindMarkers(hires,`ident.1` = ident1,`ident.2` = ident2,only.pos = TRUE)
    RNAmarkers <- rownames(diff_SCT %>% arrange(p_val_adj) %>% head(100))
    
    plotData <- genes %>% mutate(pos = floor(((start+end)/2+resolution/2)/resolution)*resolution) %>% right_join(cpg_diff) %>% mutate(log_p_val= -log10(p_val))
    #plotData <- genes %>% mutate(pos = floor(((start+end)/2)/resolution)*resolution) %>% right_join(cpg_diff) %>% mutate(log_p_val= -log10(p_val))
    plotData <- plotData %>% full_join(tibble(name=RNAmarkers,type="Marker")) %>% filter(!is.na(chr))
    plotData$type[is.na(plotData$type)] <- "Other"
    
    return(plotData)
}
plotCpGDiffVocanoPlot <- function(plotData, plotTitle,text_x_Pos=0.04,name_show_num =15,scale_factor = 1){
    plotData <- plotData %>% mutate(mean_diff = scale_factor * mean_diff)
    plot <- ggplot(plotData %>% filter(type=="Other")) + geom_point(aes(x=mean_diff,y=log_p_val),color='grey',alpha=0.1) + 
    geom_point(data = plotData %>% filter(type=="Marker"),aes(x=mean_diff,y=log_p_val),shape=1) + 
    ggtitle(plotTitle)  + theme(plot.title = element_text(hjust = 0.5))  + 
    geom_vline(aes(xintercept=mean((plotData%>% filter(type=="Marker"))$mean_diff)),linetype="longdash") + 
    geom_vline(aes(xintercept=0))  + 
    annotate("text",label = paste0("Mean diff: ",as.character( round(mean((plotData%>% filter(type=="Marker"))$mean_diff),digits = 3))), x = text_x_Pos, y = 0,size = 4 ) + 
    geom_text_repel(data = plotData %>% filter(type == "Marker") %>% arrange(desc(log_p_val)) %>% head(name_show_num),
                 aes(x=mean_diff,y=log_p_val,label=name),max.overlaps = getOption("ggrepel.max.overlaps", default = 20))

    
    print(paste0("mean difference of CpG of all gene is ",as.character( mean((plotData%>% filter(type=="Marker"))$mean_diff))))
    plot
}

# functions for analysis pairs file
# Usage:
# bintable is a dataframe with chrom-start-end as the first 3 columns

# library(furrr)
# plan(multicore, workers = 100)
# options(future.rng.onMisuse="ignore")
# intermingle_result <- future_map(colnames(E75_hires),~calcInterRatio(paste0("../../HiC_clean3/clean3_pairs/",.x,".pairs.gz"),cluster_res) %>% select(4)) %>% suppressMessages()
# intermingle_result <- bind_cols(intermingle_result)
# colnames(intermingle_result) <- str_extract(colnames(intermingle_result),"(Gas|Org)[a-zA-Z0-9]+")

calcInterRatio <- function(file,bintable){
    library(valr)
    pairs <- readPairs(file)
    pairs <- pairs %>% mutate(type = ifelse(chrom1 == chrom2 ,"intra","inter"))

    left <- pairs %>% select(chrom1,pos1,type) %>% mutate(end = pos1 + 1) %>% select(chrom1,pos1,end,type)
    names(left) <- c("chrom","start","end","type")
    right <- pairs %>% select(chrom2,pos2,type) %>% mutate(end = pos2 + 1) %>% select(chrom2,pos2,end,type)
    names(right) <- c("chrom","start","end","type")

    # betoolsr is the R wrapper of bedtools and not suitable for parallel computation, use valr instead.
    #intermingle_res <- bt.intersect(bintable %>% select(1:3) ,rbind(left,right),wa = T,wb=T,loj=TRUE) %>%
    #     group_by(V1,V2,V3) %>% summarise(count = sum(V7 == "inter") / n()) %>% ungroup()
    #names(intermingle_res) <- c("chrom","start","end","interRatio")
    intermingle_res <- bed_intersect(bintable %>% select(1:3),rbind(left,right)) %>% 
        group_by(chrom,start.x,end.x) %>% summarise(count = sum(type.y == "inter") / n())
    names(intermingle_res) <- c("chrom","start","end",file)
    intermingle_res <- bintable %>% select(1:3) %>% left_join(intermingle_res)
    return(intermingle_res)
}

# compartment strength calculation
# 1. Tanay's method and simple method

#example of load requirements
#pairsPaths <- paste0("../../HiC_clean3/clean3_pairs/",dir("../../HiC_clean3/clean3_pairs/"))
#resolution = 1000000
#read_tsv("../pileup_stage_lineage/majorgroup/processed/compartment/Emb.compartment.1m.cis.vecs.tsv") %>% mutate(type = ifelse(E1 > 0 ,"A","B")) %>% na.omit() %>% select(chrom,start,end,type) %>% write_tsv("emb.abcomp.1m.tsv",col_names = FALSE)
#defineAB <- read_table("./emb.abcomp.1m.tsv",col_names=FALSE)
#names(defineAB) <- c("chromosome","start","end","AB")

#example of use this funciton
#tanay_res <- calcCompScore(pairsPaths,defineAB,method = "tanay",threads = 200)

calcCompScore <- function(pairsPaths,defineAB,method = "tanay",threads = 100,bedtoolsPath = "~/miniconda3/envs/py3/bin/"){
    registerDoParallel(threads)
    result <- foreach(filepath=pairsPaths , .combine="rbind", .packages=c("tidyverse","bedtoolsr")) %dopar% {
        options(bedtools.path = bedtoolsPath)
        #c(filepath,calcCompartmentStrengthTanay(pairsPath = filepath,bulkAB = defineAB,resolution=1000000))
        testdata <- read_table(filepath,comment = "#",col_names = F,col_types = cols(X2 = col_character() ,X4 = col_character())) 
        if(dim(testdata)[1] < 100){
            return(NA)
        }
        testdata <- testdata %>% select(2,3,4,5) %>% mutate(pair_id = row_number())     
        names(testdata) <- c("chr1","pos1","chr2","pos2","pair_id")
        testdata <- testdata %>% mutate(pos1e = pos1 + 1,pos2e = pos2+1)  %>% filter(chr1 == chr2)
        #testData <- cbind(bt.intersect(a=testdata %>% select(chr1,pos1,pos1e),b=defineAB,wa = TRUE,wb=TRUE,loj=TRUE),
        #                 bt.intersect(a=testdata %>% select(chr2,pos2,pos2e),b=defineAB,wa = TRUE,wb=TRUE,loj=TRUE) ) %>% na.omit()
        left <-  testdata %>% select(chr1,pos1,pair_id) %>% mutate(pos11 = floor(pos1 / resolution)*resolution)
        names(left) <- c("chromosome","pos1","pair_id","start")
        left <- left %>% left_join(defineAB) %>% select(-end)

        right <-  testdata %>% select(chr1,pos2,pair_id) %>% mutate(pos21 = floor(pos2 / resolution)*resolution)
        names(right) <- c("chromosome","pos2","pair_id","start")
        right <- right %>% left_join(defineAB) %>% select(-end)
        testdata <- cbind(left,right)

        names(testdata) <- paste0("V",seq(1:10))
        testdata <- testdata #%>% filter(!(V7 - V2 < 2000000))
        testdata <- testdata %>% select(V5,V10)
        names(testdata) <- c("left","right")
        ABcount <- testdata %>% group_by(left,right) %>% summarise(count = n()) %>% filter(left != "." & right != ".") %>% arrange(left,right)
        #compartmentStrength = as.numeric(ABcount[1,3]) * as.numeric(ABcount[4,3]) / ((as.numeric(ABcount[2,3]) + as.numeric(ABcount[3,3]))**2)
        Oaa = as.numeric(ABcount[1,3])
        Oab = as.numeric(ABcount[2,3]) + as.numeric(ABcount[3,3])
        Obb = as.numeric(ABcount[4,3])
        
        if(method == "tanay"){
            T = Oaa + Oab + Obb
            Ascore = (2*Oaa + Oab)/T
            Bscore = (2*Obb + Oab)/T
            CompScore = log2(2*Ascore*Bscore*T/Oab)
            return(c(filepath,CompScore,Ascore,Bscore))
        }
        if(method == "simple"){
            CompScore = (Oaa+Obb)/(2*Oab)
            Ascore = Oaa/Oab
            Bscore = Obb/Oab
            return(c(filepath,CompScore,Ascore,Bscore))
        }     
    }
    result <- result%>% as.data.frame() %>% as_tibble()
    names(result) <- c("pairs","compScore","PA","PB")
    result$compScore <- as.numeric(result$compScore)
    result$PA <- as.numeric(result$PA)
    result$PB <- as.numeric(result$PB)
    return(result)
}

# 2. approximation by calculate gini index of scAB value
# scAB value of each cell can be provided as input:x
# return the gini index of scAB value, and A and B score

calcCompGini <- function(x, n = rep(1, length(x)), unbiased = TRUE, conf.level = NA, R = 1000, type = "bca", na.rm = FALSE){
    x <- as.numeric(x)
    x <- rep(x, n)
    if (na.rm) 
        x <- na.omit(x)
    if (any(is.na(x)) || any(x < 0)) 
        return(NA_real_)
    i.gini <- function(x, unbiased = TRUE) {
        n <- length(x)
        x <- sort(x)
        res <- 2 * sum(x * 1:n)/(n * sum(x)) - 1 - (1/n)
        if (unbiased) 
            res <- n/(n - 1) * res
        return(pmax(0, res))
    }
    if (is.na(conf.level)) {
        res <- i.gini(x, unbiased = unbiased)
    }
    else {
        boot.gini <- boot(x, function(x, d) i.gini(x[d], unbiased = unbiased), 
            R = R)
        ci <- boot.ci(boot.gini, conf = conf.level, type = type)
        res <- c(gini = boot.gini$t0, lwr.ci = ci[[4]][4], upr.ci = ci[[4]][5])
    }
    
    n=length(x)
    lcres <- Lc(x)
    Ascore = 1 - lcres$L[(length(lcres$L)/2+1):length(lcres$L)] %>% sum() / (3*length(lcres$L)/8)
    Bscore = 1 - lcres$L[1:(length(lcres$L)/2)] %>% sum() / (length(lcres$L)/8)
    return(c(res,Ascore,Bscore))
}

