require(parallel)
require(foreach)
library(data.table)
library(dplyr)
library(tidyr)
library(MASS)
source(unlist(snakemake@params[["us_func"]]))


# set input and output files
clust_p<- unlist(snakemake@input[["clust_info"]])
genome_info_p<-unlist(snakemake@input[["genome_info"]])
stb_p<-unlist(snakemake@input[["stb"]])
motupan_p<-unlist(snakemake@input[["motupan"]])
clust_filt_p<- unlist(snakemake@input[["clust_filt"]])
core_div_p<-unlist(snakemake@input[["gene_info"]])
IS<- unlist(snakemake@input[["IS"]])
tables<- paste(IS, "/output/", sep="")
file_path <- list.files(tables, pattern="_profile_SNVs.tsv.gz", full.names = T)

print(paste(".gz files: ", file_path))

outfile<- unlist(snakemake@output[["rare"]])

# Functions
format_data<- function(data_p, sk=0){
  vec_p<- unlist(strsplit(data_p, split="/"))
  file<- vec_p[length(vec_p)]
  vec_f<-  unlist(strsplit(file, split="_"))
  sample<- vec_f[1]
  
  df <- read.csv(data_p, header=T, skip = sk, sep="\t") %>% 
    dplyr::mutate("sample"=sample) %>% 
    relocate(sample)
  
  return(df)
}

find_genome<- function(row){
  gene<- unlist(strsplit(row[2], split=";"))
  
  df<- data.frame("gene"= gene) 
  
  return(df)
}

format_motupan_file<- function(file, sk){
  print(file)
  df<- read.csv(file, skip=sk, header=T, sep="\t") %>%
    filter(type=="core", mean_copy_per_genome==1) %>%
    dplyr::select(genomes, genes)
  
  formatted_df<- do.call(rbind, apply(df, 1, find_genome))
  
  return(formatted_df)
}


# open cluster info
print("opening clustering info")
clust<- read.csv(clust_p, header=T, sep="\t")
print(head(clust))
print("\n")

#open genome info
print("opening genome info")
genome_info<-  read.csv(genome_info_p, header=T, sep="\t") %>%
  dplyr::mutate("id"= get_mtdata(mtdata_path=clust_p,
                          "id", corr="genome", genome))%>%
  relocate(species, .after = genome) %>%
  relocate(genus, .after = species) %>%
  relocate(id, .after = genus)

genome_info_filt <- genome_info %>% filter(coverage_median>=5 & breadth >= 0.5) %>%
  group_by(sample) %>%
  dplyr::mutate("frequency"=filtered_read_pair_count/sum(filtered_read_pair_count)) %>%
  relocate(frequency, .before = coverage) %>%
  group_by(sample, genus) %>%
  dplyr::mutate(freq_inG=coverage/sum(coverage))

## get a list of id and the samples where they are found
sam_id <- genome_info_filt %>% dplyr::mutate(sam_ids=paste(sample, id, sep="-")) %>% pull(sam_ids)
sam_id <- unique(sam_id)

print(head(genome_info_filt))
print("\n")

# open core genes info
print("opening core genome info")
stb<- read.csv(stb_p, header=F, sep="\t")
colnames(stb)<- c("scaffold", "genome")

all_genes_df<- fread(core_div_p) %>% left_join(., stb, by="scaffold") %>%
  dplyr::mutate("species"=get_mtdata(mtdata_path=clust_p, "species", corr="genome", genome),
         "genus"= get_mtdata(mtdata_path=clust_p, "genus", corr="genome", genome),
         "id"=get_mtdata(mtdata_path=clust_p, "id", corr="genome", genome),
         "secondary_cluster"=get_mtdata(mtdata_path=clust_p, "secondary_cluster", corr="genome", genome),
         "secondary_cluster_n"=get_mtdata(mtdata_path=clust_p, "secondary_cluster_n", corr="genome", genome))

gene_to_genome<- setNames(all_genes_df$genome, all_genes_df$gene)
gene_to_genome<- gene_to_genome[unique(names(gene_to_genome))]

# format motupan data
print(motupan_p)
motupan_df<-do.call(rbind, lapply(motupan_p,function(x) format_motupan_file(x, sk=16)))%>% 
  dplyr::mutate(genome=gene_to_genome[gene]) %>% drop_na()
rownames(motupan_df)<- 1:nrow(motupan_df)

# Get a list of core genes according to mOTUpan
core_genes_motupan<- unique(motupan_df$gene)

# Update the gene file by adding info on core and accessory plus filtering for ids that are found in the sample 
all_genes_df<- all_genes_df %>% filter(genus %in% get_core()) %>%
  dplyr::mutate(type_motupan=ifelse(gene %in% core_genes_motupan, "core", "accessory"),
         sam_ids=paste(sample, id, sep="-"))

all_genes_df_filt<- all_genes_df %>%
  filter(sam_ids %in% sam_id, type_motupan=="core")

print(head(all_genes_df_filt))
print("\n")

# perform rarefaction
print("performing rarefaction...")

# Set file path
add<- all_genes_df_filt %>% dplyr::select(gene, id, gene_length)

# Initialize empty vector to store the number of unique SNVs per sample
snv_counts_df <- data.frame(sample = numeric(), id = character(), n_snvs = numeric())
ges<- unique(all_genes_df_filt$gene)

# Loop through all files
get_snvs_rarefaction<- function(gz_file, gene_wanted=ges, all_snvs){
  file<- basename(gz_file)
  sample<- unlist(strsplit(file, split = "_"))[1]
  print(paste("counting unique SNVs in sample:", sample))
  
  # Read in file
  sample_data <- fread(gz_file, header = T, sep= "\t", select = c("gene", "position","position_coverage", "class")) %>% 
    filter(gene %in% gene_wanted, position_coverage>=5, grepl("SNV", class)) %>%
    left_join(.,add, by="gene") %>%
    dplyr::mutate(sample=sample, SNV_position=paste(gene, position, sep="_")) %>%
    dplyr::select(sample, SNV_position, gene, gene_length, id)%>%
    relocate(sample) %>% distinct()
  
  # Group the data by bacteria_species and count the number of unique SNVs in each group
  #sample_data[,.(gene, gene_length, n_snvs = uniqueN(setdiff(SNV_position, unique(all_snvs)))), by = id]
  
  snv_counts_species <- sample_data %>% 
    filter(SNV_position %in% setdiff(SNV_position, all_snvs))%>%
    group_by(id)%>%
    dplyr::summarise(n_snvs=uniqueN(SNV_position),
                     sample=unique(sample),
                     len_genes=sum(gene_length),
                     perc_snv=n_snvs/len_genes)
  
  # Add all the SNV_position to the all_snvs vector
  all_snvs <- unique(c(all_snvs, sample_data$SNV_position))
  snv_counts_species$total_snvs<- length(all_snvs)
  print(paste("new SNVs=", length(all_snvs)))
  rm(sample_data)
  
  # Append the number of unique SNVs for each species in the sample to the snv_counts vector
  final_df <- snv_counts_species
  
  return(list(final_df, all_snvs))
}


#create the cluster
n.cores <- parallel::detectCores()
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK",
  outfile=""
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)


snv_counts_df<- NULL
set.seed(3033244)
x<- foreach(i=1:10)%dopar%{
  suppressPackageStartupMessages(library(foreach))
  # Initialize empty vector to store all the SNV_position across all samples
  all_snvs_l <- c()
  print(paste("\n#########rarefaction bootstrap", i))
  
  #shuffle the sample list
  file_path<- sample(file_path, length(file_path))
  
  y<- foreach(f = file_path)%dopar%{
    suppressPackageStartupMessages({library(dplyr)
    library(parallel)
    library(data.table)
    library(tidyr)
    library(MASS)})
    
    a<- get_snvs_rarefaction(f, all_snvs=all_snvs_l)[[1]]
    all_snvs_l<- get_snvs_rarefaction(f, all_snvs=all_snvs_l)[[2]]
    a$group<- i
    snv_counts_df<- rbind(snv_counts_df, a)
  }
  
  
}

snv_counts_df<- bind_rows(x) %>% distinct()
print("saving rarefaction file")
write.table(snv_counts_df, outfile, quote=F, row.names = F, sep="\t")
print("DONE!")

