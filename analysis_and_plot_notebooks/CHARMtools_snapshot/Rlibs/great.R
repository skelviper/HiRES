# wrappers for GREAT analysis

library(rGREAT)

runGREAT <- function(bed,species = "mm10",bg = NULL,method = "simple"){
    bed <- bed %>% select(1:3)
    names(bed) <- c("chr","start","end")
    job = submitGreatJob(bed,bg = bg,species = species,request_interval=1,version = "4.0.4",
        rule="basalPlusExt")#, adv_upstream = 1, adv_downstream = 1, adv_span = 1)
    tb = getEnrichmentTables(job,download_by = 'tsv')

    return(tb)
}