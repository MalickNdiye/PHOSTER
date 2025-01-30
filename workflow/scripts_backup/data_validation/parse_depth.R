library(tidyverse)
library(vroom)

if("data.table" %in% rownames(installed.packages())){
  library(data.table)
} else{
  install.packages("data.table", repos = "http://cran.us.r-project.org")
}
## TODO include the first output for cov profiles

#filelist<- unlist(snakemake@input[["pb_depth"]])
restats_files<- unlist(snakemake@input[["refstats"]])

summarise_df<- function(filename){
  bsname=basename(filename)
  sample=strsplit(bsname, "_")[[1]][1]


  df<-vroom(filename, col_names = F, delim="\t")
  head(df)
  colnames(df)<- c("contig", "locus", "depth")
  df_sum<- df %>% group_by(contig)%>%
    summarise("IQR_depth"=IQR(depth), "median_depth"=median(depth),
              "norm_depth"=IQR(depth)/median(depth), "max_depth"=max(depth)) %>%
    mutate("sample"=sample) %>%
    mutate("sam_type"=ifelse(grepl("P", sample), "virome", "bacteria"),
  "hit_type"=ifelse(grepl("VC_", contig), "GB_VCs", "Bacteria")) %>%
    filter(max_depth>0) %>%
    relocate(sample, .before = contig)

  return(df_sum)
}

open_refstats<- function(filename){
  bsname=basename(filename)
  sample=strsplit(bsname, "_")[[1]][1]

  df<-  read.table(filename, header=F, sep="\t")
  colnames(df)<- c("name", "perc_unambiguousReads", "unambiguousMB",	"perc_ambiguousReads",	"ambiguousMB",
                   "unambiguousReads",	"ambiguousReads",	"assignedReads", "assignedBases")

  df_new<- df %>% mutate("sample"=sample) %>%
    mutate("sam_type"=ifelse(grepl("P", sample), "virome", "bacteria"),
           "hit_type"=ifelse(grepl("VC_", name), "GB_VCs", ifelse(grepl("Amel", name),"Honeybee", "Bacteria"))) %>%
             relocate(sample, .before = name)

           return(df_new)
}

# print(filelist[1])
# final_df<- summarise_df(filelist[1])
# for(i in 2:length(filelist)){
#   print(filelist[i])
#   df<- summarise_df(filelist[i])
#   final_df<- rbind(final_df, df)
# }

#final_df<- do.call("rbind", lapply(filelist, summarise_df)) #WARNING very slow
refstats_df<-do.call("rbind", lapply(restats_files, open_refstats))

#write_delim(final_df, snakemake@output[[1]], delim="\t")
write_delim(refstats_df, snakemake@output[[1]], delim="\t")
