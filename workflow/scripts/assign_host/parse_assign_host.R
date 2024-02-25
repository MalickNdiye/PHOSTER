#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)
library(data.table)
source("scripts/useful_func.R")

option_list = list(
  make_option(c("-m", "--spacers_metadata"), type="character", default=NULL, 
              help="metadata on spacers", metavar="character"),
  make_option(c("-d", "--drep_data"), type="character", default=NULL, 
              help="dereplication data", metavar="character"),
  make_option(c("-b", "--blastout"), type="character", default=NULL, 
              help="spacers blast result", metavar="character"),
  make_option(c("-p", "--prophages"), type="character", default=NULL, 
              help="fastANI output", metavar="character"),
  make_option(c("-c", "--bact"), type="character", default=NULL, 
              help="bacteria metadata", metavar="character"),
  make_option(c("-v", "--virus_metadata"), type="character", default=NULL, 
              help="vOTU assignation", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="output phage-host file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



# set paths
spacer_mtdata_p<- opt$spacers_metadata
virus_mtadata_p<- opt$virus_metadata
clust_info_p<- opt$bact
blastout_p<- opt$blastout
votu_p<- opt$drep_data
fastani_p<-opt$prophages

# open Virus mtdata to get quality
v_mtadata<- fread(virus_mtadata_p)
virus_quality<- setNames(v_mtadata$checkv_quality, v_mtadata$contig_id)
quality_thresh<- c("Medium-quality", "High-quality", "Complete")

# open spacers metadata
print("open and format spacers metadata")

spacers_mtdata<- read.csv(spacer_mtdata_p, sep="\t", header = T) %>%
  ungroup() %>%
  mutate(SPACER_ID=paste(Strain, Orientation, ShortID, sep="_"),
         spacer_len=as.integer(nchar(SpacerSeq)))

spacers_to_retain<- spacers_mtdata$SPACER_ID
sp_len<- setNames(spacers_mtdata$spacer_len, spacers_mtdata$SPACER_ID)

# open balst results
print("open and format blast results")
get_bact<- function(x){
  bact<- unlist(strsplit(x, split="_\\+_|_-_|_ND_"))[1]
  return(bact)
}

blastout<- fread(blastout_p) %>%
  mutate(spacer_len=sp_len[SPACER_ID],
         bacteria=unlist(lapply(SPACER_ID, get_bact))) %>%
  filter(SPACER_ID %in% spacers_to_retain) %>%
  mutate(genus=get_mtdata(mtdata_path=clust_info_p,
                          "genus", corr="genome", bacteria),
         id=get_mtdata(mtdata_path=clust_info_p,
                       "id", corr="genome", bacteria),
         true_mm=spacer_len - alignement_length + mismatch) %>%
  filter(abs(true_mm)<3)

phage_host_spacers<- blastout %>% select(Query, bacteria, id, genus) %>%
  rename("virus"=Query) %>%
  mutate(detection="CRISPR",
         state="unceratin") %>%
  distinct(virus, bacteria, .keep_all = T)

# open fastani results
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
         state=ifelse(breadth>0.8, "integrated", "uncertain")) %>% 
  filter(breadth>=0.5)

phage_host_prophages<-fastani %>% select(Query, bacteria, id, genus, state) %>%
  mutate(detection="Homology") %>%
  rename("virus"=Query) %>%
  distinct(virus, bacteria, .keep_all = TRUE)

# combine data
print("combine data")
phage_host<- rbind(phage_host_spacers, phage_host_prophages)%>%
  group_by(virus, bacteria) %>%
  mutate(
    detection=ifelse(all(c("CRISPR", "Homology") %in% detection), "Both", unique(detection)),
    state=ifelse("integrated" %in% state, "integrated", "uncertain")
         ) %>%
  distinct(virus, bacteria, .keep_all = T)

# add de-replication data
votu=read.csv(votu_p, sep="\t", header=T) %>% select(genome, representative, vOTU)
votu_remaining<- votu %>% filter(!genome %in% phage_host$virus) %>%
  mutate(bacteria="unclassified",
         genus="unclassified",
         id="unclassified",
         state=NA,
         detection=NA)%>%
  rename("virus"=genome)

phage_host_afinal<- left_join(phage_host, votu, by=c("virus"="genome")) 
phage_host_final<- rbind(phage_host_afinal, votu_remaining)
phage_host_final$checkv_quality<- virus_quality[phage_host_final$virus]

ph<- phage_host_final %>%
  filter(checkv_quality %in% quality_thresh)%>%
  group_by(vOTU) %>%
  reframe(host=length(detection[!is.na(detection)]))

n_dect<- length(ph$vOTU[ph$host>0])
n_undect<- length(ph$vOTU[ph$host==0])

frac_det<-round((n_dect/(n_dect+n_undect))*100)

print(paste0(frac_det, "% of viruses have been linked to at least one host in the database"))

# save results
## check if directory of outfile exists, if not create it
if(!dir.exists(dirname(opt$outfile))){
  dir.create(dirname(opt$outfile))
}

write.table(phage_host_final, file=opt$outfile, sep="\t", row.names = F, quote = F)
