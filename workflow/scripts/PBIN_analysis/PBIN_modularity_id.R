### get libraries
library(lpbrim)
library("optparse")
library(tidyverse)
library(data.table)
source("scripts/useful_func.R")

# get arguments
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="input PBIN", metavar="character"),
  make_option(c("-r", "--genomemat"), type="character", default=NULL, 
              help="PBIN genome level", metavar="character"),
  make_option(c("-b", "--bact"), type="character", default=NULL, 
              help="bacterial metadtata", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="prefix", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

real_mat<- opt$infile # real PBIN matrix

output_dir<- opt$outdir
dir.create(output_dir, showWarnings = FALSE)

get_mods<- function(mat, prefix, levs){
  ### This function uses the LP_BRIM algorithm to detect modules in a PBIN matrix ###
  
  # set_outfiles
  mod_out<- file.path(output_dir,paste(prefix, "mod_obj.RData", sep="_"))
  bact_out<- file.path(output_dir,paste(prefix, "bacterial_modules.txt", sep="_"))
  vir_out<- file.path(output_dir,paste(prefix, "viral_modules.txt", sep="_"))
  
  # open matrix file
  mod_mat<- read.csv(mat, sep="\t", header=T, row.names = 1) %>%
    as.matrix()
  print(colnames(mod_mat))
  print(rownames(mod_mat))
  
  # make sure it is a binary matrix
  mod_mat[mod_mat>0]<- 1
  
  # find modules
  mod_obj<-findModules(mod_mat, sparse=F, iter=200)
  save(mod_obj, file=mod_out)
  
  # create dataframes of modules assigment for bacteria and viruses
  mods_bact<- as.data.frame(mod_obj$S) %>% 
    rownames_to_column("id") %>%
    filter(!grepl("vOTU|vBin", id))%>%
    mutate(genus= get_mtdata(mtdata_path=opt$bact,
                            "genus", corr="id", id))%>%
    group_by(id, genus) %>%
    gather(module, val, -id, -genus)%>%
    filter(val == 1) %>%
    filter(id!="Bombilactobacillus mellifer cl.-49_1") %>% # has only one phage-host pair
    group_by(id) %>%
    mutate(int_module= names(sort(table(module), decreasing = TRUE)[1]))%>%
    ungroup()%>%
    mutate(int_module=  with(., match(int_module, unique(int_module))),
          module=as.character(module),
          int_module=as.character(int_module)) %>%
    arrange(as.numeric(int_module))%>%
    group_by(int_module)%>%
    mutate(module_name=paste(unique(genus), collapse = "_"),
          module_name=paste(module_name, int_module, sep=" "))

    genomes<- fread(opt$genomemat) %>% pull("V1")

    mods_bact<-fread(opt$bact) %>%
      select(genome, id) %>%
      filter(genome %in% genomes) %>%
      left_join(., mods_bact, by="id")
    
    new_mod<- mods_bact %>% group_by(module) %>%
      count(int_module) %>%
      top_n(1)
    
  new_mod<- setNames(new_mod$int_module,as.character(new_mod$module))
  names(new_mod)
  
  mods_vir<- as.data.frame(mod_obj$S) %>% 
    rownames_to_column("genome") %>%
    filter(grepl("vOTU|vBin", genome)) %>%
    group_by(genome) %>%
    gather(module, val, -genome) %>%
    filter(val == 1) %>%
    dplyr::select(-val) %>%
    mutate(module=as.character(module),
          int_module=new_mod[module]) %>%
    arrange(as.numeric(int_module))%>%
    arrange(as.numeric(int_module))
  
  # Save dataframes
  write.table(mods_vir, vir_out, quote=F, sep="\t")
  write.table(mods_bact, bact_out, quote=F, sep="\t")
  
}

# Set levels of the core members
genus_phyl=c("Gilliamella", "Frischella", "Snodgrassella", "Commensalibacter", "Bartonella" ,
             "Bifidobacterium", "Lactobacillus", "Bombilactobacillus")

set.seed(50)
get_mods(real_mat, as.character(opt$p), genus_phyl)
