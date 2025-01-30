library(tidyverse)

if("data.table" %in% rownames(installed.packages())){
  library(data.table)
} else{
  install.packages("data.table", repos = "http://cran.us.r-project.org")
}

file_list<- unlist(snakemake@input)
file_list<-paste(file_list, ".clstr", sep="")

summarise_clst<- function(filename){
  clstr <- read.csv(filename, sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)

  clstr2 <- clstr
  n = nrow(clstr)

  x = 0
  numbers_only <- function(x) !grepl("\\D", x)
  for (row in c(1:n)) {
    if (numbers_only(clstr2[row,1]) == TRUE) {
      clstr2[row,1] <- x}
    else {NULL}
    x <- clstr2[row,1]
  }


  clstr.sums <- data.frame(dplyr::count(clstr2,V1))

  switch <- clstr.sums[1,2]
  clstr3 <- cbind(clstr2[1], clstr)

  clstr4 <- clstr2[-which(clstr2$V2 == ""), ]
  clstr4


  strip<- function(str, pos, split){
    a<- unlist(strsplit(str, split))
    return(a[pos])
  }

  clstr5<- clstr4 %>% mutate(V1=unlist(lapply(V1, strip, pos=2, split=">")),
                             length=unlist(lapply(V2, strip, pos=1, split=",")),
                             length=unlist(lapply(length, strip, pos=1, split="nt")),
                             tmp=unlist(lapply(V2, strip, pos=2, split=">")),
                             contig=unlist(lapply(tmp, strip, pos=1, split="... ")),
                             identity=unlist(lapply(tmp, strip, pos=2, split="... ")),
                             direction=unlist(lapply(identity, strip, pos=1, split="/")),
                             direction=ifelse(direction=="*", "+", unlist(lapply(direction, strip, pos=2, split="at "))),
                             identity=ifelse(identity=="*", "100", unlist(lapply(identity, strip, pos=2, split="/"))),
                             identity=unlist(lapply(identity, strip, pos=1, split="%")),
                             sample=unlist(lapply(tmp, strip, pos=1, "_")),
                             kmer_cov=unlist(lapply(tmp, strip, pos=7, "_")),
                             kmer_cov=unlist(lapply(kmer_cov, strip, pos=1, "... "))) %>%
    select(-c(tmp, V2)) %>%
    rename(cluster=V1) %>%
    relocate(sample, .before=cluster) %>%
    group_by(cluster) %>%
    mutate(cluster_size=length(cluster))

  return(clstr5)
}

a<- summarise_clst(file_list[1])
b<- summarise_clst(file_list[2])

write_delim(a, unlist(snakemake@output[1]), delim="\t")
write_delim(b, unlist(snakemake@output[2]), delim="\t")
