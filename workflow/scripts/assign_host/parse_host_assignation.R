library(tidyverse)
library(data.table)
source("./scripts/useful_func.R")

# set path to data
print("set path to data")
print(paste("path to data:", unlist(snakemake@input["spacers_mtdata"])))
print(paste("path to data:", unlist(snakemake@input["clust_filtered"])))
print(paste("path to data:", unlist(snakemake@input["all_blastout"])))
print(paste("path to data:", unlist(snakemake@input["binning_data"])))
print(paste("path to data:", unlist(snakemake@input["pro"])))

spacer_mtdata_p<- paste(unlist(snakemake@input["spacers_mtdata"]), "PerSpacer_CRISPR_v1.csv", sep="/")
clust_info_p<- unlist(snakemake@input["clust_filtered"])
blastout_p<- unlist(snakemake@input["all_blastout"])
binning_data_p<-unlist(snakemake@input["binning_data"])
fastani_p<- unlist(snakemake@input["pro"])

# read and format binning data
print("read and format binning data")
binning_data<- read.csv(binning_data_p, header=T, sep="\t") %>%
  filter(filtered=="yes")
bin_cont<- setNames(binning_data$bin, binning_data$contig)
bin_cont<-bin_cont[unique(names(bin_cont))]

# read and format spacer metadata
print("read and format spacer metadata")
get_spacer_len<- function(seq){
  len<-length(unlist(strsplit(seq, split="")))
  return(len)
}


spacer_mtdata<- read.csv(spacer_mtdata_p, header=T, sep=";") %>% 
  filter(EvidenceLvl==4) %>%
  mutate(spacer_length=unlist(lapply(SpacerSeq, get_spacer_len)),
         secondary_cluster=
           get_mtdata(mtdata_path=clust_info_p,
                      "secondary_cluster", corr="genome", Strain),
         secondary_cluster_n=
           get_mtdata(mtdata_path=clust_info_p,
                      "secondary_cluster_n", corr="genome", Strain),
         SPACER_ID=paste(Strain, Orientation, ShortID, sep="_"),
         genus=get_mtdata(mtdata_path=clust_info_p,
                          "genus", corr="genome", Strain),
         id=get_mtdata(mtdata_path=clust_info_p,
                       "id", corr="genome", Strain),
         genome_type=ifelse(grepl("MAG", Strain), "MAG", "Isolate")) %>%
  distinct(SpacerSeq, id, .keep_all = TRUE)

# read and format blastout
print("read and format blastout")
blastout<- fread(blastout_p, header=T) %>%
  left_join(., spacer_mtdata, by="SPACER_ID") 
blastout[blastout==""]<- NA

## filter blastout for match with max 2 true mistmatches 
print("filter blastout for match with max 2 true mistmatches")
blastout_filt<- blastout %>% 
  mutate(spacer_len=ifelse(!is.na(SPACER_LENGTH),
                           SPACER_LENGTH,spacer_length),
         Strain=ifelse(!is.na(Strain), Strain,GENEBANK_ID),
         SpacerSeq=ifelse(!is.na(SpacerSeq), SpacerSeq,SPACER),
         id=ifelse(!is.na(SPECIES), SPECIES,id),
         SpacerNb=ifelse(!is.na(POSITION_INSIDE_LOCUS), POSITION_INSIDE_LOCUS,SpacerNb),
         genus=ifelse(!is.na(GENUS), GENUS,genus),
         true_mm=spacer_len - alignement_length + mismatch,
         genome_type=ifelse(is.na(genome_type),
                            "CRISPRdb", genome_type),
         spacer_origin=ifelse(genome_type=="CRISPRdb", "CRISPRdb",
                              "mine")) %>%
    filter(true_mm < 3, gap==0) %>%
  select(- c(GENEBANK_ID, SPECIES, ORGANISM_NAME, GENUS,
             FAMILY, ORDER,SPACER_LENGTH, ShortID,SPACER, COUNT_SPACER, POSITION_INSIDE_LOCUS,
             true_num_mismatch, spacer_length)) %>%
  relocate(id, genus, Strain, Query, Hit_nr, identity, true_mm)
rm(blastout)

## get only spacers from my genomes
print("get only spacers from my genomes")
blastout_filt_mine<- blastout_filt %>%
  filter((spacer_origin=="mine" | !genus %in% get_core()),
         genus!="Unknown")

## spacers host is the full table, also with contig that have not be classified
print("spacers host is the full table, also with contig that have not be classified")   
bin_red<- binning_data %>% select(contig, quality, bin, length)

spacers_host<- blastout_filt_mine %>% 
  select(Query, id, genus, Strain,  genome_type, spacer_origin) %>%
  full_join(., bin_red, by=c("Query"="contig")) %>%
  mutate(binned=ifelse(grepl("bin", bin), "yes", "no"),
         Strain=ifelse(is.na(Strain), "unclassified", Strain),
         genus=ifelse(is.na(genus), "unclassified", genus),
         id=ifelse(is.na(id), "unclassified", id)) %>%
  group_by(bin) %>%
  mutate(bin_size=length(unique(Query))) 


# open and format fastani host assignation
print("open and format fastani host assignation")
format_name_fastani<- function(name){
  name_list <- unlist(strsplit(name, split="/"))
  n<- length(name_list)
  new_name<- name_list[n]
  
  if(grepl("fasta", new_name)){
    new_name<- unlist(strsplit(new_name, split=".fasta"))[1]
      
  }
  
  return(new_name)
}


fastani<- fread(fastani_p, sep="\t", 
                   col.names =c("Query", "bacteria", "ANI",
                                "frag_n", "total_query_frag")) %>%
  filter(ANI>=90) %>%
  mutate(Query=unlist(lapply(Query, format_name_fastani)),
         bacteria=unlist(lapply(bacteria, format_name_fastani)),
         breadth=frag_n/total_query_frag,
         id=get_mtdata(mtdata_path=clust_info_p,
                       "id", corr="genome", bacteria),
         genus=get_mtdata(mtdata_path=clust_info_p,
                          "genus", corr="genome", bacteria),
         state=ifelse(breadth>0.8, "integrated", "uncertain")) 

# Merge Host Data
print("Merge Host Data")
spacers_host_red<- spacers_host %>% ungroup() %>%
  select(Query, Strain, genus, id, bin) %>%
  mutate(assign_type=ifelse(genus=="unclassified", "None", "CRISPR")) %>%
  rename("virus"=Query, "bacteria"=Strain)

fastani_red<-fastani %>%
  select(Query, bacteria, genus, id)%>%
  mutate(assign_type="identity",
         bin=bin_cont[Query]) %>%
  rename("virus"=Query)

phage_host<- rbind(spacers_host_red, fastani_red) %>% group_by(virus, bacteria, genus, id) %>%
  summarise(detection=ifelse(length(unique(assign_type)[!assign_type=="None"])<=1, unique(assign_type), "Both"),
            bin=bin_cont[virus])

det<- phage_host %>% ungroup() %>% group_by(virus, genus) %>%
  summarise(score=ifelse(genus=="unclassified", 0, 1)) %>% group_by(virus) %>%
  summarise(sum_s=sum(score),
            det=ifelse(sum_s==0, "no", "yes"))
detected<- 100-length(det$det[det=="no"])/length(det$det)*100

print(paste(detected, "% of contigs have been assigned to a host"))

# save files
print("save files")
write.table(spacer_mtdata, unlist(snakemake@output["formatted_spacers"]), sep="\t", quote = F, row.names = F )

write.table(blastout_filt, unlist(snakemake@output["blastout_filt_complete"]), sep="\t", quote = F, row.names = F )

write.table(blastout_filt_mine, unlist(snakemake@output["blastout_formatted"]), sep="\t", quote = F, row.names = F )

write.table(spacers_host, unlist(snakemake@output["spacers_host"]), sep="\t", quote = F, row.names = F )

write.table(fastani, unlist(snakemake@output["prophages_host"]), sep="\t", quote = F, row.names = F )

write.table(phage_host, unlist(snakemake@output["phage_host"]), sep="\t", quote = F, row.names = F )
