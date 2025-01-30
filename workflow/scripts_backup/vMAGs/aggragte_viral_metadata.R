#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)
library(data.table)

option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="viral identification metadata", metavar="character"),
  #make_option(c("-b", "--binning_data"), type="character", default=NULL, 
   #           help="binning data", metavar="character"),
  make_option(c("-t", "--trimming_data"), type="character", default=NULL,
              help="prophage_trimming", metavar="character"),
  make_option(c("-q", "--quality_data"), type="character", default=NULL,
              help="quality filtering data", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="output summary file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# set table path
## These file have been created in order of appearence, but the name of the 
## contig changes. So, the point of this script is to reconcile all the info
## to the final name of the contigs.
mtdata_p<- opt$metadata
trimming_p<- opt$trimming_data
quality_p<- opt$quality_data

# select sample
sam<- unlist(strsplit(basename(quality_p), split = "_"))[1] 


mtdata<- fread(mtdata_p) %>% 
  filter(sample==sam) %>%
  dplyr::select(new_contig, type, prophage)
binning<-  data.frame(contig=character(), bin_name=character())

mtadata_bin<- left_join(mtdata, binning, by=c("new_contig"="contig"), multiple="all") %>%
  mutate(bin_name=ifelse(is.na(bin_name), new_contig, bin_name),
         type=ifelse(prophage=="yes", "lysogenic", type)) %>%
  group_by(bin_name) %>%
  reframe(type=ifelse("lysogenic" %in% type, "lysogenic", unique(type)),
          prophage=ifelse("yes" %in% prophage, "Yes", "No"))

trimming<- fread(trimming_p) %>%
  dplyr::select(contig, new_name, provirus) 


if(nrow(trimming)==0){
  trimming<- data.frame(contig=character(), new_name=character(), provirus=character())
}

metadata_trim<- full_join(mtadata_bin, trimming, by=c("bin_name"="contig")) %>%
  mutate(new_name=ifelse(is.na(new_name), bin_name, new_name))%>%
  select(-bin_name) %>%
  rename("contig"=new_name) %>%
  group_by(contig)%>%
  mutate(provirus=ifelse(is.na(provirus), "No", provirus),
         prophage=ifelse(provirus=="Yes", "Yes", prophage)) %>%
  mutate(type=ifelse(prophage=="Yes", "lysogenic", type)) %>%
  select(- provirus)


quality<- fread(quality_p) 
metadata_final<- full_join(quality, metadata_trim, by=c("contig_id"="contig")) %>%
  select(contig_id, contig_length, prophage, type, gene_count,
         viral_genes, host_genes, checkv_quality, completeness_method, retained, reason) %>%
  mutate(type=ifelse(is.na(type),"lytic", type),
  sample=sam) %>%
  relocate(sample)


# write outfile
write.table(metadata_final, opt$outfile, sep="\t", quote = F, row.names = F)