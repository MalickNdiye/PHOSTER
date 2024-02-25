#!/usr/bin/env Rscript
library("optparse")
library(tidyverse)

option_list = list(
  make_option(c("-i", "--trim"), type="character", default=NULL, 
              help="CheckV contamination folder", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default="out.txt", 
              help="output file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$trim)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# open & format data
## checkV contamination
trim_p<- paste(file.path(opt$trim, "contamination.tsv"))
trim<- read.csv(trim_p, sep="\t", header=T)

trim_red<- trim %>% 
  select(contig_id, contig_length, provirus, region_types, region_lengths, region_coords_bp) 

# Functions
split_coord<- get_entries<- function(str){
  # this function split the coordinates in the checkV file
  coords<- unlist(strsplit(unname(str), split = "-"))
  
  return(coords)
}

format_row<- function(x){
  # This function format the checkV data to have one file per contigs once the contig is decontaminated
  if(x["provirus"]=="No"){
    df<- data.frame("contig"=x[1], "contig_length"=x[2],"provirus"=x[3], 
                    "start"=NA, 
                    "end"=NA,
                    "prophage_length"=NA)
    return(df)
  }
  
  regions<- unlist(strsplit(unname(x["region_types"]), split = ","))
  viral_reg<- which(regions=="viral")
  n<- length(viral_reg)
  
  viral_length<- unlist(strsplit(unname(x["region_lengths"]), split = ","))[viral_reg]
  viral_coord<- unlist(strsplit(unname(x["region_coords_bp"]), split = ","))[viral_reg]
  
  df<- data.frame("contig"=vector(), "contig_length"=vector(),
                  "provirus"=vector(), "start"=vector(), 
                  "end"=vector(),
                  "prophage_length"=vector())
  for(i in 1:n){
    print(regions)
    new_df<- data.frame("contig"=x["contig_id"], "contig_length"=x["contig_length"],
                        "provirus"=x["provirus"], "start"=split_coord(viral_coord[i])[1], 
                        "end"=split_coord(viral_coord[i])[2],
                        "prophage_length"=viral_length[i])
   
    df<- rbind(df, new_df)
  }
  
  return(df)
}

# format trimming data
trim_format<- do.call(rbind, apply(trim_red, 1, format_row))
rownames(trim_format)<- 1:nrow(trim_format)

trim_final<- trim_format %>%
  group_by(contig) %>%
  mutate(n=n(),
         new_name= ifelse(n>1, paste0(contig, "_", row_number()), contig)) %>%
  select(-n)

# select only proviruses
checkv_proviruses<- trim_final %>% filter(provirus=="Yes")

# save file
write.table(checkv_proviruses, opt$outfile, quote = F, row.names = F, sep="\t")


