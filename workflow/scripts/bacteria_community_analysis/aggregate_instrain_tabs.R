library(tidyselect)
source(unlist(snakemake@params[["us_func"]]))


###
# This script takes all the data tables in the output directories of instrain for every samples and concatenat the results.
###

# get path clust info
clust_info_p<- unlist(snakemake@input[["tax"]])

# get instrain profiles data tables
profile_list<- unlist(snakemake@input[["IS"]])
data_tab_addon<- "/output/"

output_list=paste(profile_list, data_tab_addon, sep="")

all_tabs<- unlist(lapply(output_list, list.files))

# Format data data tables
format_data<- function(data_p, sk=0){
  file<- basename(data_p)
  vec_file<- unlist(strsplit(file, split="_"))
  sample<- vec_file[1]
  
  df <- read.csv(data_p, header=T, skip = sk, sep="\t") %>% 
    dplyr::mutate("sample"=sample) %>% 
    relocate(sample)
  
  return(df)
}

format_data_list<- function(path_list, pattern, get_taxa_info=F){
  files_list<- path_list[grepl(pattern, path_list)]
  
  if(get_taxa_info){
    df<- do.call(rbind, lapply(genome_info_l, format_data))%>%
      mutate("species"= get_mtdata(mtdata_path="clust_filtered_winners.tsv", "species", corr="genome", genome),
             "genus"=get_mtdata(mtdata_path="clust_filtered_winners.tsv", "genus", corr="genome", genome),
             "SDP"= get_mtdata(mtdata_path="clust_filtered_winners.tsv", "SDP", corr="genome", genome),
             "sec_clust"= get_mtdata(mtdata_path="clust_filtered_winners.tsv", "secondary_cluster", corr="genome", genome))
  }else{
    df<- do.call(rbind, lapply(genome_info_l, format_data))
  }
  
  return(df)
}

# concatenating instrain outputs
print("--- Concatenating instrain outputs...")
genome_info<- format_data_list(all_tabs, "genome_info.tsv", get_taxa_info=T)
gene_info<- format_data_list(all_tabs, "gene_info.tsv")
mapping_info<- format_data_list(all_tabs, "mapping_info.tsv")
scaffold_info<- format_data_list(all_tabs, "scaffold_info.tsv")

print("--- Saving concatenated tables...")
write.table(genome_info, file=unlist(snakemake@output[["genome_info"]]), quote=F, row.names = F, sep="\t")
write.table(gene_info, file=unlist(snakemake@output[["gene_info"]]), quote=F, row.names = F, sep="\t")
write.table(mapping_info, file=unlist(snakemake@output[["mapping_info"]]), quote=F, row.names = F, sep="\t")
write.table(scaffold_info, file=unlist(snakemake@output[["scaffold_info"]]), quote=F, row.names = F, sep="\t")

print("--- DONE!")

