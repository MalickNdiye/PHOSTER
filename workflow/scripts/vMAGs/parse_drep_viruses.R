#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)
library(data.table)

option_list = list(
  make_option(c("-i", "--drep_data"), type="character", default=NULL, 
              help="dRep_directory", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="output summary file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# open file
data_p<-file.path(opt$drep_data, "data_tables/Cdb.csv")
drep_data<-read.csv(data_p, header=T, sep=",") %>%
  select(genome, greedy_representative, primary_cluster, secondary_cluster) %>%
  mutate(vOTU=paste0("vOTU_", as.numeric(factor(secondary_cluster)))) %>%
  rename("representative"=greedy_representative) %>%
  mutate(representative=gsub("False", "FALSE", representative),
         representative=gsub("True", "TRUE", representative),
         genome=gsub(".fasta", "", genome))

print(paste0("Number of genomes: ", nrow(drep_data)))
print(paste0("Number of vOTU: ", nrow(drep_data[drep_data$representative=="TRUE",])))

# write file
write.table(drep_data, opt$outfile, sep="\t", row.names=F, quote=F)