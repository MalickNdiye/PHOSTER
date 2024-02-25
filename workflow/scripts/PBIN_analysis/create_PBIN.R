#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)
library(data.table)
source("scripts/useful_func.R")

option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="input phage_host file", metavar="character"),
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

# Open phage_host data
phage_host<- read.csv(opt$infile, sep="\t", header = T)

# create PBIN matrix
get_int_df <- function(df, quality_thersh, genus_phyl, target_genus) {
  
  interaction_df<- df %>% filter(genus %in% target_genus,
                                 checkv_quality %in% quality_thersh) %>%
    mutate(genus=as.factor(genus),
           det_value=ifelse(
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
write.table(blastout_matrix_iso, file=file.path(opt$outdir, "blastout_matrix_iso.tsv"), quote=F, sep="\t")


