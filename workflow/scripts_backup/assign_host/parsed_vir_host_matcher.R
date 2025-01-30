library("optparse")
library(tidyverse)
library(data.table)
library(reshape2)
source("scripts/useful_func.R")

option_list = list(
  make_option(c("-i", "--vhm"), type="character", default=NULL, 
              help="list of vhm directories", metavar="character"),
  make_option(c("-c", "--clust_info"), type="character", default=NULL, 
              help="data on bacteria phylogeny", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

indirs=unlist(strsplit(opt$vhm, split=","))

format_vhm<- function(vhm_dir){
  path<- file.path(vhm_dir, "d2star_k6.csv")
  print(paste("reading", path, "..."))

    vhm<- fread(path) %>% melt() %>% filter(!is.na(value))  %>% select(-"V700")
    print(head(vhm))
    colnames(vhm)<- c("virus", "bacteria", "d2star")
    vhm<- vhm %>% 
    mutate(d2star=round(d2star, 2)) %>%
    filter(d2star<0.25)

    return(vhm)

}

vhm<- do.call(rbind, lapply(indirs, format_vhm))
vhm<- vhm %>% distinct(virus, bacteria, .keep_all = TRUE) %>%
    mutate(detection="Homology",
            state="uncertain",
            bacteria=gsub(".fa|.fna", "", bacteria),
            virus=gsub(".fasta", "", virus),
            id=get_mtdata(mtdata_path=opt$clust_info,
                       "id", corr="Bin.Id", bacteria),
            genus=get_mtdata(mtdata_path=opt$clust_info,
                          "genus", corr="Bin.Id", bacteria),)%>% 
            #dplyr::select(-d2star) %>%
            relocate(virus, bacteria, id, genus, detection, state)

# save output
write.table(vhm, opt$output, sep="\t", quote=F, row.names=F)