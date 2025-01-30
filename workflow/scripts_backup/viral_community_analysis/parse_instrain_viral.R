### get libraries
library("optparse")
library(tidyverse)
library(data.table)
source("scripts/useful_func.R")

# get arguments
iqrs<- snakemake@input[["iqr"]]
indirs<- snakemake@input[["IS"]]
drep<- snakemake@input[["drep"]]
metadata<- snakemake@input[["metadata"]]
outdir<- snakemake@output[["outdir"]]

# quality contigs
quality_thresh=c("Medium-quality", "High-quality", "Complete")

# viral metadata vectors
viral_mtdata<- fread(metadata) 
head(viral_mtdata)
lifestyle_phages<- setNames(viral_mtdata$type, viral_mtdata$contig_id) 
qual_phages<- setNames(viral_mtdata$checkv_quality, viral_mtdata$contig_id) 
gen_len<- setNames(viral_mtdata$contig_length,viral_mtdata$contig_id)
gen_len<- gen_len[!duplicated(names(gen_len))]

# drep data
drep_data<- fread(drep)
head(drep_data)

#vOTUs to retain
vOTU_qual_keep<- drep_data %>% filter(representative==T, checkv_quality %in% quality_thresh) %>%
  pull(vOTU)

cluster_phage<- setNames(drep_data$vOTU, drep_data$genome)
cluster_phage<- cluster_phage[!duplicated(names(cluster_phage))]

# format iqr data
format_iqr<- function(file){
    sample<-unlist(strsplit(basename(file), split="_"))[1]
    df<- fread(file) 
    colnames(df)<- c("genome", "iqr")
    df<- df %>% mutate(sample=sample)

    return(df)
}

iqr_data<- do.call(rbind, lapply(iqrs, format_iqr))

print(head(iqr_data))


# open genome data
# append "output" to each file in opt$indir 
directories<- unlist(lapply(indirs, function(x) file.path(x, "output")))
vircom_info_list<- unlist(lapply(directories, function(x) list.files(x, pattern="genome_info.tsv", full.names=TRUE)))

format_genome_info<- function(file){
  sample<-unlist(strsplit(basename(file), split="_"))[1]
  
  df<- fread(file) %>%
    mutate(sample=sample,
           genome=gsub(".fasta", "", genome)) %>%
    relocate(sample)
  
  return(df)
}
print(vircom_info_list)

vircom_info<- do.call(rbind, lapply(vircom_info_list, format_genome_info) ) %>% as.data.frame() %>%
    left_join(., as.data.frame(iqr_data), by=c("sample", "genome")) %>%
    mutate("quality"= qual_phages[genome],
         type=lifestyle_phages[genome],
         cluster_id=cluster_phage[genome],
         IQRm=iqr/coverage_median,
         IQRm=ifelse(IQRm==Inf, 0, IQRm)) %>%
  group_by(sample)%>%
  mutate(total_reads=sum(filtered_read_pair_count)) %>% 
  ungroup() %>%
  mutate(min_reads=min(total_reads)) %>%
  group_by(sample)%>%
  mutate(norm_fact=min_reads/total_reads,
         rarefied_reads=round(norm_fact*filtered_read_pair_count,0)) %>%
  relocate(cluster_id, quality, .after=genome) 

mean(vircom_info$IQRm)

vircom_info_filt_qual<- vircom_info %>% filter(breadth>=0.70, cluster_id %in% vOTU_qual_keep)%>%
  group_by(sample) %>%
  mutate("actual_length"=length*breadth,
         norm_reads=filtered_read_pair_count/actual_length,
         "frequency"=norm_reads/sum(norm_reads),
         genome_length=gen_len[genome]) %>%
  relocate(sample, genome, cluster_id, quality, frequency, coverage, coverage_median, IQRm)

# save vircom_info_filt_qual
# if directory does not exist, create it
if(!dir.exists(outdir)){
  dir.create(outdir)
}
write.table(vircom_info, file.path(outdir, "vircom_info.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(vircom_info_filt_qual, file.path(outdir, "vircom_info_filt_qual.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
