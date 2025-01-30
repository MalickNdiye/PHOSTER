#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)

option_list = list(
  make_option(c("-i", "--checkv"), type="character", default=NULL, 
              help="CheckV quality file", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$checkv)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


# set variables
quality_ret<- c("Low-quality", "Medium-quality", "High-quality", "Complete")
min_len<-10000

# open checkv data
checkv_p<- opt$checkv
checkv<-read.csv(checkv_p, sep="\t", header=T) %>%
  mutate(retained=1,
         reason="retained")

# Filtering criteria
## 1. length>10k and (contig  shorter are very low quality)
## 2. Low- to Complete-quality (low quality are used only for mapping but excluded from analyses)
## 3. kmer frequency <1.1 (too high indicates contamination)
checkv_filt<- checkv %>%
  mutate(retained=ifelse(contig_length<min_len, 0, retained),
         reason=ifelse(contig_length<min_len, "length", reason)) %>%
  mutate(retained=ifelse(checkv_quality %in% quality_ret, retained, 0),
         reason=ifelse(checkv_quality %in% quality_ret,  reason, "quality")) %>%
  mutate(retained=ifelse(kmer_freq<1.1, retained, 0),
         reason=ifelse(kmer_freq<1.1,  reason, "contamination"))



retained_n<- length(checkv_filt$contig_id[checkv_filt$retained==T])
all_contigs<- length(checkv_filt$contig_id)
perc_ret<- round((retained_n/all_contigs)*100)

print(paste0(retained_n, " contigs were retained as likely phages out of the initial ", all_contigs, "--> ", perc_ret, "% retained contigs."))

# write output
# check if output directory exists, if not create it
if (!dir.exists(dirname(opt$outfile))){
  dir.create(dirname(opt$outfile))}
write.table(checkv_filt, opt$outfile, sep = "\t", row.names = F,  quote=F)