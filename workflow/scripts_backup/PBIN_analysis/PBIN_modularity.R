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
  set.seed(569)
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
  mods_bact<- as.data.frame(mod_obj$S) %>% rownames_to_column("genome") %>%
    filter(!grepl("vOTU|vBin", genome))%>%
    mutate(genus= get_mtdata(mtdata_path=opt$bact,
                             "genus", corr="genome", genome),
           id=get_mtdata(mtdata_path=opt$bact,
                         "id", corr="genome", genome))%>%
    mutate(genus=ifelse(is.na(genus), get_mtdata(mtdata_path=opt$bact,
                                                      "genus", corr="Bin.Id", genome), genus),
           id=ifelse(is.na(id), get_mtdata(mtdata_path=opt$bact,
                                                      "id", corr="Bin.Id", genome), id)) %>%
    group_by(genome, genus, id) %>%
    gather(module, val, -id, -genome, -genus) %>%
    filter(val == 1) %>%
    select(-val) %>% 
    mutate(module=ifelse(id=="Gilliamella sp. cl.-30_1", unique(.$module[.$id=="Gilliamella apicola cl.-35_1"]), module)) %>% # Put all the Gilliamella sp. cl.-30_1 in the same module as the other Gillimellas because all MAGs of this species are basically identical, i.e. the module is an artifact
    mutate(module=ifelse(id=="Bifidobacterium coryneforme  cl.-40_1", unique(.$module[.$id=="Bifidobacterium polysaccharolyticum cl.-41_1"]), module)) %>% # Put all the Bifidobacterium coryneforme  cl.-40_1 in the same module as the other bifido because it creates a module for itslef even though it has most phages shared with other bifidos
    filter(id!="Bombilactobacillus mellifer cl.-49_1") %>% #  has only one phage-host pair, not treatable
    arrange(factor(genus, levels=levs), id, module) %>%
    group_by(id) %>%
    mutate(int_module= names(sort(table(module), decreasing = TRUE)[1])) %>%
    ungroup()%>%
    mutate(int_module=  with(., match(int_module, unique(int_module))),
           module=as.character(module),
           int_module=as.character(int_module)) %>%
    arrange(as.numeric(int_module))%>%
    group_by(int_module)%>%
    mutate(module_name=paste(unique(genus), collapse = "_"),
           module_name=paste(module_name, int_module, sep=" "))
  
  # count genomes in moodules
  new_mod<- mods_bact %>% group_by(module) %>%
    count(int_module) %>%
    top_n(1)
  
  new_mod<- setNames(new_mod$int_module,as.character(new_mod$module))
  names(new_mod)
  
  mods_vir<- as.data.frame(mod_obj$S) %>% rownames_to_column("genome") %>%
    filter(grepl("vOTU|vBin", genome)) %>%
    group_by(genome) %>%
    gather(module, val, -genome) %>%
    filter(val == 1) %>%
    select(-val) %>%
    mutate(module=as.character(module),
           int_module=new_mod[module]) %>%
    ungroup()%>%
    mutate(int_module=ifelse(is.na(int_module),1, int_module), # add the phages that belong to Gilliamella sp. cl.-30_1 to the module gilliamella 1
      int_module=ifelse(genome=="vOTU_863" & prefix=="allGenomes", 7, int_module),
      int_module=ifelse(genome=="vOTU_863" & prefix=="Isolates", 6, int_module))%>%  #add the only phage that belong to Bifido coryneforme to the module of the other bifidos (a bit rough)
    arrange(as.numeric(int_module))%>%
    arrange(as.numeric(int_module))
  
  # Save dataframes
  write.table(mods_vir, vir_out, quote=F, sep="\t")
  write.table(mods_bact, bact_out, quote=F, sep="\t")
  
}

# Set levels of the core members
genus_phyl=c("Gilliamella", "Frischella", "Snodgrassella", "Commensalibacter", "Bartonella" ,
             "Bifidobacterium", "Lactobacillus", "Bombilactobacillus")

get_mods(real_mat, as.character(opt$p), genus_phyl)
