#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)
library(data.table)
source("scripts/useful_func.R")

option_list = list(
  make_option(c("-m", "--spacers_metadata"), type="character", default=NULL, 
              help="metadata on spacers", metavar="character"),
  make_option(c("-d", "--DR"), type="character", default=NULL, 
              help="direct repeat metadata", metavar="character"),
  make_option(c("-c", "--bact"), type="character", default=NULL, 
              help="bacteria metadata", metavar="character"),
  make_option(c("-o", "--outspacer"), type="character", default=NULL,
              help="output spacer file", metavar="character"),
  make_option(c("-e", "--outdr"), type="character", default=NULL,
              help="output DR file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# bacteria metadata
clust_info_p<- opt$bact
# open spacers metadata
spacers_mtdata<-read.csv(opt$spacers_metadata, sep=";", header=T)  %>% dplyr::filter(EvidenceLvl==4)

array_qual<- setNames(spacers_mtdata$EvidenceLvl, spacers_mtdata$ArrayID)

# open DR cluster
print("Filtering DR metadata...")
DR_clust_raw<- read.csv(opt$DR, sep="\t") %>% select(-representative)
DR_clust<- DR_clust_raw %>%
  rowwise()%>%
  mutate(ArrayID=unlist(strsplit(contig_id, split=".fna_|.fa_"))[2],
         genome=unlist(strsplit(contig_id, split=".fna_|.fa_"))[1]) %>%
  ungroup() %>%
  mutate(id=get_mtdata(mtdata_path=clust_info_p,
                                      "id", corr="Bin.Id", genome),
         genus=get_mtdata(mtdata_path=clust_info_p,
                          "genus", corr="Bin.Id", genome),
         EvidenceLvl=array_qual[ArrayID],
         type=ifelse(grepl("DR1", contig_id), "metaG", "genome"),
         genome_type=ifelse(grepl(".fna", contig_id), "isolate", "MAG"),
         genome_type=ifelse(grepl("DR1", contig_id), "metaG", genome_type)) %>%
  filter(!(type=="genome" & is.na(EvidenceLvl)))%>% # remove low confidence arrays
  group_by(cluster_id)%>%
  mutate(n_genus_iso=length(unique(genus[!is.na(genus) & genome_type=="isolate"])),
         n_genus=length(unique(genus[!is.na(genus)])),
         cluster_size=length(cluster_size)) %>%
  select(-sample) %>%
  relocate(genome, genus, id) %>% ungroup()

# filter spacers mtadata
print("Filtering spacers metadata...")
spacers_mtdata_filt<- spacers_mtdata %>%  dplyr::filter(EvidenceLvl==4) %>%
  mutate(id=get_mtdata(mtdata_path=clust_info_p,
                       "id", corr="genome", Strain),
         genus=get_mtdata(mtdata_path=clust_info_p,
                          "genus", corr="genome", Strain))  %>%
  relocate(genus, id, .after=Strain)%>% ungroup()

sort(unique(spacers_mtdata_filt$id))

# save file
## check if folders exist, if not create them
if(!dir.exists(opt$outspacer)){
  dir.create(opt$outspacer)
}
if(!dir.exists(opt$outdr)){
  dir.create(opt$outdr)
}

write.table(spacers_mtdata_filt, file.path(opt$outspacer,"spacers_metadata_filtered.tsv"), sep="\t", row.names = F, quote = F)
write.table(DR_clust, file.path(opt$outdr, "DR_metadata_filtered.tsv"), sep="\t", row.names = F, quote = F)
