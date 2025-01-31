#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)
library(data.table)
source("scripts/useful_func.R")

option_list = list(
    make_option(c("-i", "--indir"), type="character", default=NULL, 
              help="MAGS directory", metavar="character"),
    make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="MAGs Metadata", metavar="character"),
    make_option(c("-b", "--bactcom"), type="character", default=NULL,
                help="bacteria_community data", metavar="character"),
    make_option(c("-o", "--outdir"), type="character", default="out.txt", 
              help="output dir  all core", metavar="character"),
     make_option(c("-g", "--outgen"), type="character", default="out.txt", 
              help="output dir genera", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# open bacterial community data
print("Opening bacterial community data")
bactcom_p<- file.path(opt$bactcom, "all_genome_info_filtered.tsv")
bactcom<- fread(bactcom_p)

## filter genomes form the core microbota that are present in at least one sample
bactom_present_core<- bactcom %>% filter(genus %in% get_core())
id_core_present<- bactom_present_core %>% pull(unique(id))

# open MAGs metadata
print("Opening MAGs metadata")
metadata<- fread(opt$metadata)
metadata_core_present<- metadata %>% filter(id %in% id_core_present)
genome_core_present<- metadata_core_present %>% pull(unique(genome))
genome_core_present_list<- split(metadata_core_present$genome, metadata_core_present$genus)

# copy files in genome_core_present_list to outdir/core_present
outdir<- opt$outdir
if(!dir.exists(file.path(outdir))){
  dir.create(file.path(outdir), recursive = TRUE)
}
# copy files
lapply(genome_core_present_list, function(x) file.copy(from=file.path(opt$indir, x), to=file.path(outdir, "core_present")))


# copy files in genome_core_present_list to outdir/"name of the list"
outdir<- opt$outgen
n<- names(genome_core_present_list)
lapply(n, function(x) {
  if(!dir.exists(file.path(outdir, x))){
    dir.create(file.path(outdir, x), recursive = TRUE)
  }
  # copy files
  lapply(genome_core_present_list[[x]], function(y) file.copy(from=file.path(opt$indir, y), to=file.path(outdir, x)))
}
)