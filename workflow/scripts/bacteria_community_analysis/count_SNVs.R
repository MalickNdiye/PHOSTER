library(data.table)
library(tidyverse)
library(optparse)

# Parse arguments
option_list = list(
  make_option(c("-i", "--instrain"), type="character", default=NULL,
              help="Input SNV file"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output SNV file"),
  make_option(c("-t", "--genome_info"), type="character", default=NULL,
              help="Genome info file"),
  make_option(c("-s", "--scaffold_to_genome"), type="character", default=NULL,
              help="Scaffold to genome file"),
  make_option(c("-d", "--drep"), type="character", default=NULL,
              help="drep_genome info")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# add .gz to the input file, if the file is compressed otherwise keep the name as is
# to uderstand if the file is compressed, check i the file + .gz exists
if(file.exists(paste0(opt$instrain, ".gz"))){
  opt$instrain<- paste0(opt$instrain, ".gz")
}

# open genome info
genome_length<- fread(opt$drep) %>%
  dplyr::select(genome, length)

# get sample name from input file
sample_name<- strsplit(basename(opt$instrain),"_")[[1]][1]

writeLines(paste("Processing sample", sample_name))

# Open genome Info file
genome_info<- fread(opt$genome_info) %>%
  dplyr::select(genome, id)

# open scaffold-to-genome file
stb<- fread(opt$scaffold_to_genome, header = F, col.names = c("scaffold", "genome"))

# Open SNV file
snvs<- fread(opt$instrain) %>%
  filter(position_coverage >= 5, class!="SNS") %>% # get SNVs with at least 5x coverage
  left_join(., stb, by="scaffold") %>%
  left_join(., genome_info, by="genome") %>%
  left_join(., genome_length, by="genome") %>%
  mutate(SNV_ID=paste(scaffold, position, sep="_"),
        sample=sample_name)%>%
  relocate(sample, genome, id, length, SNV_ID)
head(snvs)

# save file as compressed tsv
writeLines(paste("\nWriting output file", opt$output))
gzfile<- gzfile(opt$output, "w")
write.table(snvs, gzfile, sep="\t", quote=F, row.names=F, col.names=T)
close(gzfile)

writeLines(paste("Finished processing sample", sample_name))
 
 





