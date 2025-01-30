#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)
library(data.table)
source("scripts/useful_func.R")

option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="input phage_host file", metavar="character"),
  make_option(c("-b", "--bact_com"), type="character", default=NULL, 
              help="bacteria community data filtered", metavar="character"),
  make_option(c("-d", "--dRep"), type="character", default=NULL, 
              help="bacteria community data filtered", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# set global variables
genus_phyl=c("Gilliamella", "Frischella", "Snodgrassella", "Commensalibacter", "Bombella", "Bartonella" ,
             "Bifidobacterium", "Lactobacillus", "Bombilactobacillus")
target_genus<- get_core()
quality_thersh<- c("Medium-quality", "High-quality", "Complete")

# Get vOTUs with at least one member satifiying quality condition
dRep<- fread(opt$dRep)
vOTU_qual_keep<- dRep %>% filter(representative==T, checkv_quality %in% quality_thersh) %>%
  pull(vOTU)
vOTU_qual_keep<- unique(vOTU_qual_keep)

# Get bacteria that are actually detected in the samples
bact_list<- read.csv(file.path(opt$bact_com, "all_genome_info_filtered.tsv"), sep="\t", header = T) %>%
  pull(unique(id))
  print(bact_list)

# Open phage_host data
phage_host<- read.csv(opt$infile, sep="\t", header = T)

# create PBIN matrix
get_int_df <- function(df, quality_thersh, genus_phyl, target_genus) {
  
  interaction_df<- df %>% filter(genus %in% target_genus,
                                 vOTU %in% vOTU_qual_keep,
                                 id %in% bact_list) %>%
    mutate(genus=as.factor(genus),
           det_value=ifelse( # det alue of 1 for CRISPR, 2 for Homology and 3 for Both
             grepl("Both", detection), 3,
             ifelse(grepl("Homology", detection), 2,
                    ifelse(grepl("CRISPR", detection), 1, 0)))) %>%
    group_by(vOTU, bacteria, genus, id)%>%
    summarise(det_value=ifelse(sum(unique(det_value))>3,
                               3,
                               sum(unique(det_value))))%>%
    group_by(vOTU) %>%
    mutate(host_range_id=length(unique(id)),
           host_range_bact=length(unique(bacteria)))%>%
    group_by(bacteria) %>%
    mutate(predating_phages=length(unique(vOTU)))%>%
    group_by(id) %>%
    mutate(total_phages_clst=length(unique(vOTU)),
           mean_phages_clst=mean(predating_phages))%>%
    ungroup()%>%
    arrange(factor(genus, levels=genus_phyl), mean_phages_clst, predating_phages,  desc( host_range_bact))
  
  return(interaction_df)
  
}

interaction_df<- get_int_df(phage_host, quality_thersh, genus_phyl, target_genus)

blastout_matrix<- interaction_df %>%
  pivot_wider(names_from = vOTU, id_cols = bacteria, values_from=det_value, values_fill=0) %>% 
  column_to_rownames("bacteria") %>%
  as.matrix()

blastout_matrix_id<- interaction_df %>%
  distinct(id, vOTU, .keep_all = T) %>%
  mutate(det_value=1)%>%
  pivot_wider(names_from = vOTU, id_cols = id, values_from=det_value, values_fill=0) %>% 
  column_to_rownames("id") %>%
  as.matrix()

interaction_df_iso<- get_int_df(phage_host, quality_thersh, genus_phyl, target_genus) %>% filter(grepl(".fna", bacteria))  %>%  ungroup()

blastout_matrix_iso<- interaction_df_iso %>% 
  pivot_wider(names_from = vOTU, id_cols = bacteria, values_from=det_value, values_fill=0) %>% 
  column_to_rownames("bacteria") %>%
  as.matrix()

# check weather outidr exists and save both matrices
if(!dir.exists(opt$outdir)){
  dir.create(opt$outdir)
}

write.table(blastout_matrix, file=file.path(opt$outdir, "blastout_matrix.tsv"), quote=F, sep="\t")
write.table(blastout_matrix_id, file=file.path(opt$outdir, "blastout_matrix_id.tsv"), quote=F, sep="\t")
write.table(blastout_matrix_iso, file=file.path(opt$outdir, "blastout_matrix_iso.tsv"), quote=F, sep="\t")


