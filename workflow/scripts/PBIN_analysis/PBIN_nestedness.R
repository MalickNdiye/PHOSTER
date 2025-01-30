### get libraries
library(tidyverse)
library(data.table)
library(vegan)
library("optparse")
source("scripts/useful_func.R")
source("scripts/PBIN_analysis/bipartite_network_analysis.R")

# get arguments
option_list = list(
    make_option(c("-p", "--PBIN"), type="character", default=NULL, 
              help="PBINs directory", metavar="character"),
    make_option(c("-m", "--modules"), type="character", default=NULL, 
              help="modules directory", metavar="character"),
    make_option(c("-i", "--iter"), type="character", default=NULL,
                help="number of iterations", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



mat_p<- file.path(opt$PBIN, "blastout_matrix.tsv")
bact_mods_p<- file.path(opt$mod, "allGenomes_bacterial_modules.txt")
vir_mods_p<- file.path(opt$mod,"allGenomes_viral_modules.txt")


# main_mat
main_mat<- fread(mat_p)%>% 
column_to_rownames("V1") %>%
as.matrix()
main_mat[main_mat>1]<- 1

# get modules dataframes
bact_mods<- fread(bact_mods_p) 

mods_id<- setNames(bact_mods$module_name, bact_mods$id)
mod_names<- setNames(bact_mods$module_name, bact_mods$int_module)
mod_names<- mod_names[!duplicated(names(mod_names))]

# Viral IMs
vir_mods<- fread(vir_mods_p) %>%
  mutate(module_name=mod_names[int_module]) 

vir_mods_vec<- setNames(vir_mods$module_name, vir_mods$genome)

mods_names<- unique(bact_mods$module_name)

# get module list
mod_list_bact<- split(bact_mods$genome, bact_mods$module_name)
mod_list_vir<- split(vir_mods$genome, vir_mods$module_name)
n_mods<- length(mod_list_bact)

names(mod_list_vir)
names(mod_list_bact)

# get module matrixes
mod_list<- lapply(mods_names, function(x) main_mat[mod_list_bact[[x]], mod_list_vir[[x]]])
names(mod_list)<- mods_names


# calculate nestdedness
library(parallel)
library(foreach)
library(doParallel)
# re-create test_df in parallel
cl <- makeCluster(12)
registerDoParallel(cl)
test_df <- foreach(i = 1:n_mods, .combine = rbind) %dopar% {
  library(tidyverse)
  library(data.table)
  library(vegan)
  library("optparse")
  source("scripts/useful_func.R")
  source("scripts/PBIN_analysis/bipartite_network_analysis.R")
  print("module to test:")
  print(mod_list[[i]])
  TestNest(mod_list[[i]], iter=opt$iter)
}

# if opt$out directory does not exist, create it
if (!file.exists(dirname(opt$out))) {
  dir.create(dirname(opt$out), recursive = TRUE)
}

write.csv(test_df, opt$out, row.names = FALSE, quote = FALSE, sep = "\t")

