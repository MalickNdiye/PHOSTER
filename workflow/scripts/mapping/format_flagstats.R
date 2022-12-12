# This scripts open the fed flagstats files as outputed by the command 'samtools flagstats' and returns them in a table format concatenating them
library(tidyverse)

filelist<-snakemake@input[["filelist"]]

message<- paste("----------------- concatenating the info found in the flloging files", filelist, sep=" ")
print(message)

format_flastat<- function(file){
  # Sample name of the table
  list_p<-unlist(strsplit(file, split = "/"))
  len<-length(list_p)
  sam<- list_p[len]
  sam_name<-unlist(strsplit(sam, split = "_"))[1]

  # columns name of the table
  cln<- c("sample", "total_reads", "primary", "secondary", "supplementary", "duplicates",
          "primary_duplicates", "mapped", "primary mapped",
          "paired_in_sequencing", "read1", "read2", "properly_paired",
          "both_mapped", "singletons", "mate_other_chr", "mate_other_chr_HQ")

  #open flagstas file and take relevant info
  nums<-read.csv(file, header=F, sep= " ") %>% pull(V1)
  values<-append(sam_name, nums)
  values<- values[1:(length(values)-1)]
  names(values)<- cln

  df<- data.frame(as.list(values))
}

final_df<-do.call(rbind, lapply(filelist, format_flastat))

write.table(final_df, file=snakemake@output[[1]], quote=F, sep="\t", row.names = F)
