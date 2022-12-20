# functions for cellcycle analysis
library(tidyverse)
library(patchwork)
library(ggpubr)
library(cellAlign)
library(Seurat)
library(circlize)
library(ComplexHeatmap)
library(DescTools)

cellcyclecolors = c("G0"="#762a83","G1"="#c2a5cf","Early-S"="#99d8c9","Mid-S"="#66c2a4","Late-S"="#41ae76","G2"="#238b45","M"="#ea6e34")

# functions for align two group of cell and compare cdps, adapt from cellAlign, NatMethod 2018
getCelltypeOrder <- function(obj,ct,sample_number = NULL){
    # order cellcycle using adapted Tanay's method(Nagano et al. 2017)
    # assuming cdps already save in a Seurat object slot named "cdps"
    # ct: celltype but actually is the expression, e.g. 'celltype == "mix late mesenchyme"'
    #     similarly 'stage == "E75"'
    
    mat <- obj[["cdps"]]@data %>% as.matrix() 
    if(is.null(sample_number)){
        tempOrderDF <- obj@meta.data %>% filter(eval(parse(text=ct))) 
    }
    else{
        tempOrderDF <- obj@meta.data %>% filter(eval(parse(text=ct)))%>% sample_n(sample_number)
    }
    
    #dim(tempOrderDF) %>% print()
    tempOrderDF <- tempOrderDF %>% mutate(clusterOrder = str_replace(paste0(str_replace(celltype," ","_"),"_",cluster,"_",sub_k_cluster),"_sub_",""))
    tempOrderDF$cellcycle_threshold <- factor(tempOrderDF$cellcycle_threshold,levels = c("G0","M","G1","Early-S","Mid-S","Late-S","G2"))
    tempMat <- tempOrderDF %>% select(cellname,clusterOrder)%>% left_join(mat %>% t() %>% as.data.frame() %>% rownames_to_column("cellname")) %>% 
        group_by(clusterOrder) %>% select(-cellname) %>% summarise_all(mean) %>% na.omit() %>% ungroup() %>% column_to_rownames("clusterOrder")

    tempOrderDF <- tempOrderDF %>% ungroup() %>%  mutate(varrepli = var(raw_repli_score),varnearp = var(near_p),snearp=scale(near_p),sfar=scale(farAvg))

    tempOrderDF <- tempOrderDF %>% group_by(cellcycle_threshold) %>% 
                 mutate(order = ifelse(cellcycle_threshold%in% c("Early-S"),near_p /varnearp + raw_repli_score / varrepli,
                                    ifelse(cellcycle_threshold %in% c("Mid-S","G2"), near_p /varnearp,
                                      ifelse(cellcycle_threshold %in% c("Late-S"),near_p / varnearp - raw_repli_score / varrepli,
                                            ifelse(cellcycle_threshold %in% c("G1","G0"),snearp + sfar,
                                                  ifelse(cellcycle_threshold %in% c("M","Unknown"),maxbinorder)))))) %>% arrange(cellcycle_threshold,order)
    return(tempOrderDF %>% pull(cellname))
}

align_cellcycle <- function(obj,celltype1,celltype2,numPts = 50,sample_number = NULL){
    mat <- obj[["cdps"]]@data %>% as.matrix() 
    traj_celltype1 <- getCelltypeOrder(obj = obj,ct = celltype1,sample_number = sample_number)%>% as.data.frame() %>% mutate(order = row_number() / n()) %>% deframe()
    traj_celltype2 <- getCelltypeOrder(obj = obj,ct = celltype2,sample_number = sample_number)%>% as.data.frame() %>% mutate(order = row_number() / n()) %>% deframe()

    interGlobal_celltype1 = cellAlign::interWeights(expDataBatch = mat[,traj_celltype1%>% names()], trajCond = traj_celltype1,winSz = 0.1, numPts = numPts)
    interGlobal_celltype2 = cellAlign::interWeights(expDataBatch = mat[,traj_celltype2%>% names()], trajCond = traj_celltype2,winSz = 0.1, numPts = numPts)
    interScaledGlobal_celltype1 = cellAlign::scaleInterpolate(interGlobal_celltype1)
    interScaledGlobal_celltype2 = cellAlign::scaleInterpolate(interGlobal_celltype2)

    # here both local and global alignmnet can be use, but for whole cellcycle i dont think local alignment is necessary since 
    # every celltype has G1~G2 and it's contact decay profile look similar.

    #alignment = localAlign(interScaledGlobal_celltype1$scaledData,interScaledGlobal_celltype2$scaledData,threshPercent = Thresh)
    #
    #mapping = mapRealDataLocal(alignment,intTrajQuery = interScaledGlobal_celltype1$traj, realTrajQuery = traj_celltype1,
    #                            intTrajRef = interScaledGlobal_celltype2$traj, realTrajRef = traj_celltype2)

    alignment = globalAlign(interScaledGlobal_celltype1$scaledData,interScaledGlobal_celltype2$scaledData,scores = list(query = interScaledGlobalML$traj, 
                                                     ref = interScaledGlobalEM$traj),sigCalc = F, numPerm = 20)
    mapping = mapRealDataGlobal(alignment,intTrajQuery = interScaledGlobal_celltype1$traj, realTrajQuery = traj_celltype1,
                                intTrajRef = interScaledGlobal_celltype2$traj, realTrajRef = traj_celltype2)

    #code for generate alignment plot
    ref_annotation <- mapping$refAssign %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("refAssign") %>% unnest(V1) %>% mutate(cellname = V1) %>% 
        select(-V1) %>% left_join(obj[[]] %>% select(cellname,cellcycle_threshold)) %>% group_by(refAssign)%>% 
        mutate(cellcycle_ref = cellcycle_threshold[which.max(n())],metaNodeRef=refAssign) %>% ungroup() %>% select(metaNodeRef,cellcycle_ref) %>% unique() 

    query_annotation <- mapping$queryAssign %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("queryAssign") %>% unnest(V1) %>% mutate(cellname = V1) %>% 
        select(-V1) %>% left_join(obj[[]] %>% select(cellname,cellcycle_threshold)) %>% group_by(queryAssign)%>% 
        mutate(cellcycle_query = cellcycle_threshold[which.max(n())],metaNodeQuery=queryAssign) %>% ungroup() %>% select(metaNodeQuery,cellcycle_query) %>% unique() 

    metaNodePt = mapping$metaNodesPt
    metaNodePt = metaNodePt[order(metaNodePt$ptQuery), ]
    metaNodePt$align = 1:nrow(metaNodePt)
    metaNodePtLong = reshape2::melt(metaNodePt[, c("ptQuery", "ptRef", "align")], id.vars = c("align"))
    metaNodePtLong = reshape2::melt(metaNodePt, id.vars = c("align", "metaNodeQuery", 
            "metaNodeRef"))

    lev2num <- function(str){
        return(which(c('G0','M','G1','Early-S','Mid-S','Late-S','G2') == str))
    }
    temp <- metaNodePtLong %>% filter(variable != "diff") %>% left_join(ref_annotation %>% group_by(cellcycle_ref) %>% mutate(rank_ref = row_number() / n())) %>%
                                                      left_join(query_annotation %>% group_by(cellcycle_query) %>% mutate(rank_query = row_number() / n())) %>% 
         rowwise() %>% mutate(shift = rank_ref + lev2num(cellcycle_ref) - rank_query - lev2num(cellcycle_query),
                              shift_old = lev2num(cellcycle_ref) - lev2num(cellcycle_query)) %>% 
        mutate(cellcycle = ifelse(variable == "ptQuery",as.character(cellcycle_query),as.character(cellcycle_ref))) %>% 
        mutate(variable = ifelse(variable == "ptQuery",str_extract(celltype1,pattern = '(?<=")[A-Za-z0-9 ]+(?=")'),str_extract(celltype2,pattern = '(?<=")[A-Za-z0-9 ]+(?=")')))

    # draw alignment plot
    temp$variable <- factor(temp$variable,levels = rev(c(str_extract(celltype1,pattern = '(?<=")[A-Za-z0-9 ]+(?=")'),str_extract(celltype2,pattern = '(?<=")[A-Za-z0-9 ]+(?=")'))))
    alignment_plot <- ggplot(temp %>% filter(variable != "diff"), aes(x = variable, y = value, group = align)) + 
            geom_line(color = "grey") + theme_bw() + geom_point(aes(color = cellcycle)) + 
            coord_flip() + ggtitle("meta-nodes alignment")  + theme_Publication()#+ scale_color_manual(values = cellcyclecolors)

    # draw histogram
    mean_num_old <- temp %>% pull(shift_old) %>% mean()
    mean_num <- temp %>% pull(shift) %>% mean()

    histogram_plot <- temp %>% gghistogram(x="shift",bins= 16) + theme_Publication() + xlab(paste0("shift (mean = ",round(mean_num,3),")")) + scale_x_continuous(limits = c(-2,2))
    #print(mean_num_old)
    #print(mean_num)
    p_return <- (alignment_plot | histogram_plot) + plot_layout(widths = c(2, 1))

    return(c(list(p_return),mean_num,list(metaNodePt),list(alignment),list(mapping)))
}

# plot cdps heatmap
plotCdpsHeatmap <- function(obj,ct){
    colors <-c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594")
    order = getCelltypeOrder(obj,ct = ct)
    tempOrderDF = obj[[]][order,]
    mat = obj[["cdps"]]@data %>% as.matrix %>% as.data.frame()
    heatmap_mat <- mat[,order] %>% as.data.frame
    heatmap_mat[heatmap_mat > 0.025] <- 0.025

    col_fun = colorRamp2(c(-1, 0, 1), c("#ffffff", "#9c76c2", "#272d74"))
    col_fun1 = colorRamp2(c(-1, 0, 1), c("#ffffff","#7eb5b4","#058786"))

    Annodf <- tempOrderDF %>% select(repli_score,G1S.Score,G2M.Score,sub_k_cluster,cluster,cellname,cellcycle_threshold,celltype)
    Annodf <- Annodf %>% group_by(cellcycle_threshold) #%>% mutate(repli_score=mean(repli_score),G1S.Score=mean(G1S.Score),G2M.Score=mean(G2M.Score))
    topAnno <- HeatmapAnnotation(df=Annodf %>% column_to_rownames("cellname")%>% select(repli_score,G1S.Score,G2M.Score,cellcycle_threshold,celltype) ,
                                col = list(repli_score = col_fun,G1S.Score = col_fun1,G2M.Score = col_fun1,cellcycle_threshold=cellcyclecolors,celltype = celltypeColors),show_annotation_name=FALSE)

    p <- heatmap_mat%>% as.matrix() %>%
        Heatmap(cluster_rows = FALSE,cluster_columns = FALSE,show_row_names = FALSE, show_column_names = FALSE,colors,
           #    column_split= tempOrderDF %>% select(cellcycle_threshold,clusterOrder),
                use_raster=TRUE,top_annotation= topAnno,column_title=str_extract(ct,pattern = '(?<=")[A-Za-z0-9 ]+(?=")'),
                column_dend_reorder = TRUE,
                   heatmap_legend_param = list(title = "contacts %"), width = 6)
    return(p)
}

# select comparable cells before use d3d
# usage:  
# sample_cell <- align_sample_cell(G0_align_res[[3]],G0_align_res[[4]],G0_align_res[[5]])
# celltype1_names_align <- sample_cell[[1]]
# celltype2_names_align <- sample_cell[[2]]
align_sample_cell <- function(metaNodePt,alignment,mapping,numPts = 50){
    localcost <- matrix(alignment$localCostMatrix,nrow = numPts) %>% as.data.frame() 
    names(localcost) <- names(mapping$refAssign)
    rownames(localcost) <- names(mapping$queryAssign)
    localcost <- localcost %>% rownames_to_column("metaNodeQuery") %>% gather(metaNodeRef,cost,-metaNodeQuery)
    sub_metanode <- metaNodePt %>% left_join(localcost) %>% group_by(metaNodeRef) %>% arrange(cost) %>% slice(1) %>% 
        group_by(metaNodeQuery) %>% arrange(cost) %>% slice(1)
    metaNodeQuery_name <- sub_metanode %>% pull(metaNodeQuery)
    metaNodeRef_name <- sub_metanode %>% pull(metaNodeRef)
    
    ref_cellname_align <- mapping$refAssign %>% t() %>% t() %>% as.data.frame() %>% rownames_to_column("metanode") %>% unnest(V1) %>% filter(metanode %in% metaNodeRef_name) %>% pull(V1) 
    query_cellname_align <- mapping$queryAssign %>% t() %>% t() %>% as.data.frame() %>% rownames_to_column("metanode") %>% unnest(V1) %>% filter(metanode %in% metaNodeQuery_name) %>% pull(V1) 
    return(list(query_cellname_align,ref_cellname_align))
}
