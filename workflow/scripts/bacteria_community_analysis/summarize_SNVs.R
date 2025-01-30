library(data.table)
library(tidyverse)
library(optparse)

# Parse arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="comma separared list of input SNV files"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output summarized table")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

writeLines(paste("Processing sample SNVs tables"))
# split input files
input_files<- unlist(strsplit(opt$input, split=","))
writeLines(paste("Input files:\n", input_files))

# rbind all input files
snvs<- do.call(rbind, lapply(input_files, fread)) 
head(snvs)

# summarize SNVs
writeLines(paste("Summarizing SNVs"))
snvs_summary<- snvs %>%
  group_by(genome) %>%
  mutate(total_snvs=uniqueN(SNV_ID)) %>%
  group_by(sample, genome) %>%
  reframe(id=unique(id),
          length=unique(length),
          total_snvs=unique(total_snvs),
          n_snvs=uniqueN(SNV_ID),
          perc_snvs=(n_snvs/length)*100,
          perc_tot_snvs=(total_snvs/length)*100) %>%
    relocate(sample, genome, id, total_snvs, n_snvs)

# save file as tsv
writeLines(paste("\nWriting output file", opt$output))
write.table(snvs_summary, opt$output, sep="\t", quote=F, row.names=F, col.names=T)
