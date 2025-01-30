#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)
library(data.table)
library(Biostrings)
source("scripts/useful_func.R")


option_list = list(
  make_option(c("-i", "--spacers_fa"), type="character", default=NULL, 
              help="FASTA of metagenomic spacers", metavar="character"),
  make_option(c("-r", "--read_count"), type="character", default=NULL, 
              help="read_count data", metavar="character"),
  make_option(c("-o", "--outspacer"), type="character", default=NULL,
              help="output spacer file", metavar="character")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# Open fasta file with Biostrings
print("---Opening Data")
fasta_p<- opt$spacers_fa
spacers<- readDNAStringSet(fasta_p)


# summarize info on spacers in a data table
print("---Summarizing Data")
break_name<- function(x, pos){
  name_list<- unlist(strsplit(x, split="_"))
  
  return(name_list[pos])
}

df<- data.frame("name"=names(spacers), len=width(spacers)) %>%
  mutate(sample=unname(sapply(name, break_name, pos=1)),
         spacer_id=unname(sapply(name, break_name, pos=2)),
         spacer_coverage=as.numeric(unname(sapply(name, break_name, pos=4))),
         sequence=unlist(as.character(spacers)),
         array=unlist(lapply(spacer_id, function(x) unlist(strsplit(x, split="SP"))[1]))) %>%
  select(-name)%>%
  relocate(sample, array, spacer_id, len, spacer_coverage, sequence)


# open read count data
print("---Opening Read Count Data")
r_count<- read.csv(opt$read_count, sep="\t", header=T) %>% select(sample, reads_postFilt) %>%
  rename(reads_postFilt="reads")

# Merge data
print("---Merging Data")
df_final<- left_join(df, r_count, by="sample") %>%
  group_by(sample)%>%
  mutate(cpm=(spacer_coverage/reads)*1000000)

# save file
print("---Saving Data")
write.table(df_final, opt$outspacer, sep="\t", row.names = F, quote = F)
