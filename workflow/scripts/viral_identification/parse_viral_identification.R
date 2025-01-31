library(tidyverse)

# Write function to format viral identification files and gives each contig a score based on the quality of the identification

format_VI<- function(path){
  path_l<- unlist(strsplit(path, split = "/"))
  path_n<- length(path_l)
  
  file=path_l[path_n]
  
  sam=ifelse(grepl("VIBRANT", file), 
             unlist(strsplit(file, split = "_"))[4],
             unlist(strsplit(file, split = "_"))[1])
  if(grepl("viralverify", path)){
    df<- read.csv(path, header=T, sep=",") %>%
      mutate(sample=sam) %>%
      relocate(sample)
  } else{
    df<- read.csv(path, header=T, sep="\t") %>%
      mutate(sample=sam) %>%
      relocate(sample)
  }
 
 
 return(df)
}

format_prophages<- function(path){
  path_l<- unlist(strsplit(path, split = "/"))
  path_n<- length(path_l)
  
  file=path_l[path_n]
  
  sam=unlist(strsplit(file, split = "_"))[1]

  df<- read.csv(path, header=T, sep="\t") %>%
    mutate(sample=sam) %>%
    relocate(sample)
  
  
  return(df)
}

# set path to viral identification results
virsorter_p<- unlist(snakemake@input["vs_score"])
vv_p<- unlist(snakemake@input["vv_score"])
vibrant_p<- unlist(snakemake@input["vib_score"])
vibrant_pro_p<- unlist(snakemake@input["vib_prophages"])
vs_pro_p<- unlist(snakemake@input["vs_prophages"])

# open all data for each identification tool and format them

## Virsorter2
virsorter<- do.call(rbind, lapply(virsorter_p, format_VI)) %>%
  rename("contig"=seqname, "len"=length, "vs_score"=max_score) %>%
  select(-"len")%>%
  filter(vs_score!="NaN") %>%
  mutate(vs_quality=ifelse(vs_score<0.5, 1, 2),
         vs_quality=ifelse((vs_score>=0.5 & vs_score<0.8), 2, vs_quality),
         vs_quality=ifelse((vs_score>=0.8), 3, vs_quality),
         ) %>%
  rowwise() %>%
  mutate(contig=unlist(strsplit(contig, split="[||]"))[1]) %>%
  distinct(sample, contig, .keep_all = T)

virsorter_pro<- do.call(rbind, lapply(vs_pro_p, format_prophages)) %>%
  filter(grepl("partial", seqname_new)) %>%
  select(sample, seqname, seqname_new, trim_bp_start, trim_bp_end) %>%
  rename("contig"=seqname, "fragment"=seqname_new,
         "start"=trim_bp_start, "stop"=trim_bp_end) %>%
  mutate(frag_len=stop-start)
  
virsorter<- left_join(virsorter, virsorter_pro, by=c("sample", "contig"))

vs_red<- virsorter %>% select(contig) %>%
  mutate(tool="Virsorter2", hit=1) %>%
  distinct()

## ViralVerify
vv<- do.call(rbind, lapply(vv_p, format_VI)) %>%
  filter(grepl("Virus|viral|too short", Prediction)) %>%
  rename("len"=Length, "vv_pred"=Prediction, "vv_score"=Score, "contig"=Contig.name) %>%
  mutate(vv_score=as.numeric(vv_score),
         vv_quality=ifelse(vv_score<5, 1,
                           ifelse(vv_score>10, 3,2)))
vv$vv_quality[is.na(vv$vv_quality)]<- 1

vv_red<- vv %>% select(contig) %>%
  mutate(tool="ViralVerify", hit=1) %>%
  distinct()

## VIBRANT
vibrant<- do.call(rbind, lapply(vibrant_p, format_VI)) %>%
  rename("vibrant_score"=Quality, "contig"=scaffold) %>%
  mutate(vibrant_quality=ifelse(vibrant_score=="low quality draft",1,
                                ifelse(vibrant_score=="medium quality draft",2,3)
                                )) %>%
  distinct(sample, contig, .keep_all=T)%>%
  rowwise() %>%
  mutate(contig=unlist(strsplit(contig, split="_fra"))[1])

vibrant_pro<- do.call(rbind, lapply(vibrant_pro_p, format_prophages)) %>%
  rename("contig"=scaffold, "fragment"=fragment,
         "start"=nucleotide.start, "stop"=nucleotide.stop,
         "frag_len"=nucleotide.length) %>%
  select(sample,contig,fragment,start,stop,frag_len)

vibrant<- left_join(vibrant, vibrant_pro, by=c("sample", "contig"))

vib_red<- vibrant %>% select(contig) %>%
  mutate(tool="VIBRANT", hit=1) %>%
  distinct()

## Merge Identified contig
VI_df<- full_join(vv, virsorter, by=c("sample", "contig")) %>%
  full_join(., vibrant, by=c("sample", "contig", "start", "stop", "fragment", "frag_len")) %>%
  rowwise()%>%
  mutate(vv_quality=ifelse(is.na(vv_quality),0, vv_quality),
         vs_quality=ifelse(is.na(vs_quality),0, vs_quality),
         vibrant_quality=ifelse(is.na(vibrant_quality),0, vibrant_quality),
         quality=mean(c(vv_quality, vs_quality, vibrant_quality), na.rm=T),
         len=as.numeric(unlist(strsplit(contig, split="_"))[4]))%>%
  ungroup()%>%
  mutate(sam_type=ifelse(grepl("P",sample), "virome", "bacteriome"),
         sam_name=factor(gsub("P|B", "", sample)),
         "prophage"= ifelse(is.na(fragment), "no", "yes")) %>%
  relocate(quality,  vv_quality, vs_quality, vibrant_quality, .after=len) %>%
  relocate(sam_type, sam_name) %>%
  select(-c(vv_pred,vv_score, Pfam.hits, dsDNAphage,
         ssDNA, vs_score, hallmark, viral, cellular)) %>%
         ungroup()

## add column named "new name" to VI_df
# the new_name will be sample_viral_contig_1, sample_viral_contig_2, etc. where the number increases within each sample
VI_df<- VI_df %>%
  group_by(sample) %>%
  mutate(new_contig=paste0(sample, "_viral_contig_", row_number())) %>%
  ungroup()%>% 
  relocate(new_contig, .after=contig) 

tool_tab<- rbind(vs_red, vv_red) %>% rbind(., vib_red)%>%
  pivot_wider(names_from="tool", values_from="hit", values_fill = 0)

# save formatted tables
write.table(VI_df,unlist(snakemake@output["all_scores"]), sep="\t", quote=F, row.names = F)
write.table(tool_tab,unlist(snakemake@output["tab"]), sep="\t", quote=F, row.names = F)
write.table(vv,unlist(snakemake@output["vv_scores"]), sep="\t", quote=F, row.names = F)
write.table(vibrant,unlist(snakemake@output["vib_scores"]), sep="\t", quote=F, row.names = F)
write.table(virsorter,unlist(snakemake@output["vs_scores"]), sep="\t", quote=F, row.names = F)
