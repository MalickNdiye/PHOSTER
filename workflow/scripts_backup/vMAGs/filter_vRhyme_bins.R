#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)

option_list = list(
  make_option(c("-d", "--bin_dir"), type="character", default=NULL, 
              help="bins directory", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="contig metadata file", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$bin_dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# get binning files
inbin=opt$bin_dir
bin_files=list.files(inbin, full.names = T)
inmember=bin_files[grepl("membership.tsv", bin_files)]
insummary=bin_files[grepl("summary.tsv", bin_files)]

# if insummary is empty, return a empty dataframe with the same columns as the output and exit with error
if (length(insummary)==0){
  bin_data_empty=data.frame(sample=character(), contig=character(), bin=character(), bin_name=numeric(),
  members=numeric(), proteins=numeric(), redundancy=numeric(), contig_length=numeric(), bin_length=numeric(),
   quality=numeric(), type=character(), prophage=character())

  if (!dir.exists(dirname(opt$outfile))){
    dir.create(dirname(opt$outfile))}
  write.table(bin_data_empty, opt$outfile, sep = "\t", row.names = F,  quote=F)
} else{
  # open and format binning tables
  bin_members=read.csv(inmember, sep="\t", header=T)
  bin_summary=read.csv(insummary, sep="\t", header=T)

  bin_data=left_join(bin_members, bin_summary, by="bin")

  # get_contig metadata
  contig_mtdata=read.csv(opt$metadata, sep="\t", header=T)

  # add info to bin_data
  bin_data_full<- left_join(bin_data, contig_mtdata, by=c("scaffold"="new_contig")) %>%
    dplyr::select(sample, scaffold, bin, members, proteins, redundancy, len, quality, type, prophage) %>%
    rename("contig_length"=len, "contig"=scaffold)%>%
    mutate(bin_name=paste0(sample, "_vBin_", bin))  %>%
    group_by(bin_name) %>%
    mutate(bin_length=sum(contig_length)) %>%
    ungroup() %>%
    relocate(bin_name, .after = bin) %>%
    relocate(bin_length, .after=contig_length)


  # filter away suspect bins
  ## remove all bin with at least one contigs where prophage=="yes"
  ## remove all bins where at least 2 contigs have type=="lysogenic", the ones with only one lysogenic contig are kept
  bin_data_filt<-bin_data_full %>%
    filter(redundancy<2)%>% 
    group_by(bin_name) %>%
    filter(!any(prophage=="yes")) %>%
    mutate(nr_lyso=sum(type=="lysogenic", na.rm = T)) %>%
    filter(nr_lyso<2) %>%
    dplyr::select(-nr_lyso) 



  # save binning file
  if (!dir.exists(dirname(opt$outfile))){
    dir.create(dirname(opt$outfile))
  }
  write.table(bin_data_filt, opt$outfile, sep = "\t", row.names = F,  quote=F)
}


