library(tidyverse)
source(unlist(snakemake@params[["us_func"]]))

## This function splits the GTDB-TK classification into diffrent taxonomical classes
get_taxa_cols<- function(entry){
  tax_l<- unlist(strsplit(entry, ";"))
  tax_l_ref<- lapply(tax_l, function(x) gsub(".*__", "", x))

  tax_l_ref[tax_l_ref==""]<- NA

  names(tax_l_ref)<- c("domain", "phylum", "class", "order", "family", "genus", "species")
  tax_l_ref$species<- gsub("__*.","", tax_l_ref$species)

  return(tax_l_ref)
}

# Set the path to data
checkm_filt_p<- unlist(snakemake@input[["checkm_filt"]])

gtdb_dir<-unlist(snakemake@input[["gtdb_dir"]])
gtdb_p<- paste(gtdb_dir, "/gtdbtk.bac120.summary.tsv", sep="")

clust_dir<-unlist(snakemake@input[["clust_dir"]])
clust_info_p<- paste(clust_dir, "/data_tables/Cdb.csv", sep="")
clust_win_p<- paste(clust_dir, "/data_tables/Wdb.csv", sep="")

mtdata_p<- unlist(snakemake@input[["mtdata"]])


# Open data
print("---Opening Data")
clust_info<-read.csv(clust_info_p, header=T)
clust_win<-read.csv(clust_win_p, header=T)
mtdata<- read.csv(mtdata_p, header=T, sep="\t")
checkm_filt<- read.csv(checkm_filt_p, header=T, sep="\t", check.names = F)
gtdb_raw<- read.csv(gtdb_p, header=T, sep = "\t")


# Format GTDB-TK data
print("---formatting GTDB-TK data")

## Here we will divide the GTDB-TK classification in different columns
gtdb<- gtdb_raw %>% group_by(user_genome) %>%
  dplyr::mutate("domain"= get_taxa_cols(classification)$domain,
                "pylum"= get_taxa_cols(classification)$pylum,
                "class"= get_taxa_cols(classification)$class,
                "order"= get_taxa_cols(classification)$order,
                "family"= get_taxa_cols(classification)$family,
                "genus"= get_taxa_cols(classification)$genus,
                "species"= get_taxa_cols(classification)$species)

## Set core genera
core=c("Bombilactobacillus" , "Commensalibacter", "Lactobacillus","Bifidobacterium","Gilliamella",  "Frischella", "Snodgrassella",  "Bartonella", "Apibacter")

## add MAGs checkm info to taxonomical classification
gtdb_QC<- left_join(checkm_filt, gtdb, by=c("Bin Id"="user_genome")) %>%
  group_by(genus) %>%
  dplyr::mutate("n_genus"= length(genus))

tax_clust<- gtdb_QC %>% select(`Bin Id`, genus, species ) %>% rename("gtdb_sp"=species)


# Format clustering data
print("---formatting clustering data")

# Format output of dRep clustering to include some taxonomical info
clust_info<- clust_info %>%
  dplyr::mutate("ref"=ifelse(genome %in% clust_win$genome, "*",""),
                "species"=get_mtdata(mtdata_p,
                                     new ="Species",
                                     df_corr = gsub(" ", "",gsub(".fna", "", genome))),
                "SDP"=get_mtdata(mtdata_p, new = "SDP",  df_corr = gsub(".fna", "", genome)),
                "Bin Id"= gsub(".fa|.fna","", genome))

# Assign a species level classification to MAGs that clusered with reference genomes or have a GtDbtk classification at the genus level
clust_assign<- clust_info %>% group_by(secondary_cluster) %>%
  dplyr::mutate(secondary_cluster_n=length(secondary_cluster),
                SDP=paste(unique(SDP[!is.na(SDP)]), collapse="/"),
                species=paste(unique(species[!is.na(species)]), collapse="/")) %>%
  left_join(., tax_clust, by="Bin Id")

clust_assign$species[(is.na(clust_assign$species) | clust_assign$species=="")]<-"unclassified"

clust_assign<- clust_assign %>%
  dplyr::mutate(species=ifelse((species=="unclassified" & is.na(gtdb_sp)), paste(genus, " sp.", sep=""),
                               ifelse((species=="unclassified" & !is.na(gtdb_sp)), gtdb_sp, species))) %>%
  select(!c("gtdb_sp"))

clust_assign$species[(is.na(clust_assign$species) | clust_assign$species=="")]<-"unclassified"


# prepare data to save
print("---prepare data to save")
gtdb_reduced<- gtdb_QC %>% ungroup()%>% select(`Bin Id`, Completeness, Contamination, `Strain heterogeneity`, `Genome size (bp)`, `N50 (contigs)`)

clust_assign_full<- left_join(clust_assign, gtdb_reduced, by="Bin Id") %>%
  dplyr::mutate("isolate"=ifelse(grepl("MAG", genome), "no", "yes"),
                genus=unlist(strsplit(species, split=" "))[1])


clust_info_full<- clust_info %>%
  dplyr::mutate("isolate"=ifelse(grepl("MAG", genome), "no", "yes"))

# # delete all cluster that are not classified at the genus level
# to_delete<- clust_assign_full %>% filter(isolate=="no") %>%
#   filter(species=="NA sp.")

clust_final<- clust_assign_full %>%  #filter(!`Bin Id` %in% to_delete$`Bin Id`) %>%
  mutate("species"=ifelse(species=="NA sp.", paste(genus, "sp.", sep=""), species),
        "species"=ifelse(species=="NAsp.", paste("Lactobacillaceae", " sp.", sep=""), species),
        "id"=paste(species, " cl.-", secondary_cluster, sep="")) %>%
  relocate("Bin Id") %>%
  relocate("species", .after=genome)

clust_final$SDP[clust_final$SDP==""]<-"unclassified"


clust_win_final<- clust_final %>% filter(ref=="*")


head(clust_final)

# save data
## This table contain the raw dRep cluster info, with some info on the taxomy
write.table(clust_info_full, unlist(snakemake@output[["clust_info"]]) , quote = F, row.names = F, sep = "\t")

## This table contains the raw dRep cluster where MAGs have been assigned a taxonomical classification in function of the references they clustered with
write.table(clust_assign_full, unlist(snakemake@output[["clust_assign"]]), quote = F, row.names = F, sep = "\t")

## This table contains the deleted dRep clusters
#write.table(to_delete, unlist(snakemake@output[["to_delete"]]), quote = F, row.names = F, sep = "\t")

## This table contains the filtered dRep cluster, with some info on the taxomy where MAGs have been assigned a taxonomical classification
write.table(clust_final, unlist(snakemake@output[["clust_final"]]), quote = F, row.names = F, sep = "\t")

## This table contains the winners of each cluster, the ones that will be used for mapping
write.table(clust_win_final, unlist(snakemake@output[["clust_win_final"]]), quote = F, row.names = F, sep = "\t")
