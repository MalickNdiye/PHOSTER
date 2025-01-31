#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)
library(data.table)
"Author: Malick Ndiaye"
option_list = list(
  make_option(c("-i", "--drep_data"), type="character", default=NULL, 
              help="dRep_directory", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="virus metadata", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="output summary file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# Open winners file
print("---Opening winners file")
winners_p<-file.path(opt$drep_data, "data_tables/Widb.csv")
winners<- fread(winners_p)

winning_genomes<- gsub(".fasta", "",winners$genome)


# get quality data
print("---Opening metadata file")
metadata_p <- opt$metadata
metdata<- read.csv(metadata_p, header=T, sep="\t") 

qual<- setNames(metdata$checkv_quality, metdata$contig_id)
qual_score<- c("Low-quality"=1, "Medium-quality"=2, "High-quality"=3, "Complete"=4)
lens<- setNames(metdata$contig_length ,  metdata$contig_id)


# open main file
print("---Opening main file")
data_p<-file.path(opt$drep_data, "data_tables/Cdb.csv")
drep_data<-read.csv(data_p, header=T, sep=",") %>%
  select(genome, primary_cluster, secondary_cluster) %>%
  mutate(vOTU=paste0("vOTU_", as.numeric(factor(secondary_cluster)))) %>%
  mutate(genome=gsub(".fasta", "", genome),
         representative=ifelse(genome %in% winning_genomes, "TRUE", "FALSE"),
         checkv_quality=qual[genome],
         qual_score=qual_score[checkv_quality],
         length=lens[genome]) %>%
  arrange(vOTU, desc(length))

print(paste0("Number of genomes: ", nrow(drep_data)))
print(paste0("Number of vOTU: ", nrow(drep_data[drep_data$representative=="TRUE",])))


# Chose best representatives
print("---Choosing best representatives")
drep_data_max <- drep_data %>%
  group_by(vOTU) %>%
  mutate(max_qual_score = max(qual_score)) %>%
  ungroup() %>%
  filter(qual_score == max_qual_score) %>% pull(genome)

vOTU_problem<- drep_data$vOTU[(!drep_data$genome %in% drep_data_max) & drep_data$representative==T]

# if there vOTU is in vOTU_problem, then change the representative to FALSE and put representative to TRUE in the best representative of the vOTU, sometimes that are multiple representatives with the same score, so chose one at random
drep_data$representative[which(drep_data$vOTU %in% vOTU_problem)]<- "FALSE"
candidates<- c()
for(j in vOTU_problem){
  drep_sub<- drep_data[drep_data$vOTU==j,]
  candidates<-append(candidates, drep_sub$genome[drep_sub$qual_score==max(drep_sub$qual_score)][1]) 
}
drep_data$representative[drep_data$genome %in% candidates]="TRUE"

drep_data<- drep_data %>% select(genome, length, checkv_quality, primary_cluster, secondary_cluster, vOTU, representative)

# write file
write.table(drep_data, opt$outfile, sep="\t", row.names=F, quote=F)