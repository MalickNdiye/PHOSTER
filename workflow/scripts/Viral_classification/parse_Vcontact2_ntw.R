# load libraries
library(tidyverse)
library(data.table)
source("./scripts/useful_func.R")

# functions
collapse_annot<- function(annot){
  annot_list<- sort(unique(annot))
  
  target_genus<- get_core()
  target_genus<- append("Apilactobaciullus", get_core())
  target_genus<- append("Bombella", get_core())
  
  if(length(intersect(target_genus, annot_list))>0){
    annot_list<- intersect(target_genus, annot_list)
  }
  
  annot_collapsed<- paste(annot_list, collapse=", ")
  annot_collapsed<- gsub(", unclassified", "", annot_collapsed)
  return(annot_collapsed)
}

check_for_core<- function(annot){
  target_genus<- get_core()
  target_genus<- append("Klebsiella", get_core())
  target_genus<- append("Apilactobaciullus", get_core())
  target_genus<- append("Bombella", get_core())
  target_genus<- append("unclassified", get_core())
  
  annot_list<- unlist(strsplit(annot, split = ", "))
  
  result<- F
  for(i in annot_list){
    if(i %in% target_genus){
      result<- T
    }
  }
  return(result)
}

# set path
phage_host_p<- unlist(snakemake@input["phage_host"])
vcont_gbg_p<- paste(unlist(snakemake@input["vcont_dir"]), "genome_by_genome_overview.csv", sep="/")
ntw_p<- paste(unlist(snakemake@input["vcont_dir"]), "c1.ntw", sep="/")

# open phage host interaction info
phage_host<- read.csv(phage_host_p, header = T, sep="\t") %>% 
  mutate(det_value=ifelse(
    grepl("Both", detection), 3,
    ifelse(grepl("identity", detection), 2,
           ifelse(grepl("CRISPR", detection), 1, 0))
  )
  )

# open and format gene-by-genome file
vcont_gbg<- fread(vcont_gbg_p) %>%
  filter(!grepl("~", Genome)) %>%
  select(Genome, "VC Status", VC, "VC Size", Quality) %>%
  rename("virus"=Genome, "VC_status" = "VC Status", "VC_size"= "VC Size", "VC_quality"=Quality) %>%
  mutate(VC=ifelse(
    VC=="",
    virus, VC),
    VC_size=ifelse(
      is.na(VC_size),
      1, VC_size
    ))

phage_host_vcont<- left_join(phage_host, vcont_gbg, by="virus")


phage_host_annot<- phage_host_vcont %>% group_by(virus, VC, VC_status, VC_size) %>%
  summarise(genus=collapse_annot(genus),
            id=collapse_annot(id),
            bacteria=collapse_annot(bacteria))

# open and format Vccontact2 Network
ntw<- fread(ntw_p, col.names = c("source_node", "end_node", "edge_weight")) %>%
  filter(!grepl("~", source_node), !grepl("~", end_node))

ntw_annot<- left_join(ntw, phage_host_annot, by=c("source_node"="virus"))

ntw_annot_filt<- ntw_annot[unlist(lapply(ntw_annot$genus, check_for_core)),]

# save files
write.table(phage_host_vcont, unlist(snakemake@output["ph_vcont"]), quote = F, sep="\t", row.names = F)

write.table(ntw_annot, unlist(snakemake@output["vcont_ntw"]), quote = F, sep="\t", row.names = F)

write.table(ntw_annot_filt, unlist(snakemake@output["vcont_ntw_filt"]), quote = F, sep="\t", row.names = F)


