################################################################################
#This parser takes in various output from the VAMB pipeline and returns some 
#formatted tables for further analyses
################################################################################
library(tidyverse)

#PATHS
## contigs mtdata
contigs_mtdata_p<- unlist(snakemake@input["contig_tab"])

## VAMB
vae_anno_p<- paste(unlist(snakemake@input["qc"]), "vambbins_aggregated_annotation.txt", sep="/")
vae_pred_p<- paste(unlist(snakemake@input["qc"]), "vambbins_RF_predictions.txt", sep="/")
dvf_pred_p<- paste(unlist(snakemake@input["all_anno"]), "all.DVF.predictions.txt", sep="/")
vamb_clst_p<- paste(unlist(snakemake@input["vae"]), "vae_clusters.tsv", sep="/")

# renamed contigs info
contigs_mtdata<- read.csv(contigs_mtdata_p, sep="\t", header=T) %>%
  filter(retained=="yes")
cont_eq<- setNames(contigs_mtdata$old_name, contigs_mtdata$new_name)


# VAMB output parsing
vamb_anno<-  read.csv(vae_anno_p, header=T, sep="\t" )
vamb_pred<-  read.csv(vae_pred_p, header=T, sep="\t" )
dvf_pred<-  read.csv(dvf_pred_p, header=T, sep="\t" ) %>%
  mutate(contig=cont_eq[name])

vamb_clst<- read.csv(vamb_clst_p, header=F, sep="\t", col.names = c("binname", "contig")) %>%
  rowwise()%>%
  mutate(contig=cont_eq[contig],
         sample=ifelse(grepl("SAMPLE", contig),
                       unlist(strsplit(contig,split="_"))[4],
                       unlist(strsplit(contig,split="_"))[1])) %>%
  left_join(., vamb_anno, by="binname")%>%
  left_join(., vamb_pred, by="binname") %>%
  left_join(., dvf_pred, by="contig") %>%
  relocate(sample, contig, name, binname) %>%
  rename("vamb_name"=name, bin_dvf_score=dvf_score)

out_vamb_p<- unlist(snakemake@output["vamb_clst"])
write.table(vamb_clst, out_vamb_p, sep="\t", quote=F, row.names=F)
# CRISPR Host-assignation


