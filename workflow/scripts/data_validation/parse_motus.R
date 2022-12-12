library(tidyverse)
library(vroom)

#
df_motus_raw  <- read.csv(snakemake@input[["motus_tab"]], sep = "\t", skip = 2, stringsAsFactors=FALSE)
df_motus_raw_count  <- read.csv(snakemake@input[["motus_count"]], sep = "\t", skip = 2, stringsAsFactors=FALSE)
colnames(df_motus_raw)[1] <- "taxonomy"
colnames(df_motus_raw_count)[1] <- "taxonomy"

df_motus_raw <- df_motus_raw %>%
              mutate(across(!c("taxonomy"), as.numeric)) %>%
              mutate(sum_ab = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>%
                filter(sum_ab > 0) %>%
                  select(!sum_ab) %>%
                    column_to_rownames("taxonomy")

df_motus_raw_count <- df_motus_raw_count %>%
              mutate(across(!c("taxonomy"), as.numeric)) %>%
              mutate(sum_count = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>%
                filter(sum_count > 0) %>%
                  select(!sum_count) %>%
                    column_to_rownames("taxonomy")


df_motus <- as.data.frame(t(df_motus_raw)) %>%
              rownames_to_column("Sample") %>%
              filter(grepl("B|ELLE", Sample)) %>%
              pivot_longer(!Sample, names_to = "motu", values_to = "rel_ab") %>%
                group_by(Sample, motu) %>%
                  mutate(Present = ifelse(rel_ab > 0, 1, 0)) %>%
                    group_by(motu) %>%
                     mutate(Prevalence_num = sum(Present),
                            Prevalence = mean(Present),
                            Sample=substr(Sample,2,1000000L)) %>%
                            group_by(Sample) %>%
                            mutate("rel_ab"=rel_ab/sum(rel_ab))

df_motus_count<- as.data.frame(t(df_motus_raw_count)) %>%
              rownames_to_column("Sample") %>%
              filter(grepl("B|ELLE", Sample)) %>%
              pivot_longer(!Sample, names_to = "motu", values_to = "read_count") %>%
                group_by(Sample, motu) %>%
                  mutate(Present = ifelse(read_count > 0, 1, 0)) %>%
                    group_by(motu) %>%
                     mutate(Prevalence_num = sum(Present),
                            Prevalence = mean(Present),
                            Sample=substr(Sample,2,1000000L))

df_combined<- left_join(df_motus, df_motus_count, by=c("Sample", "motu",
 "Present", "Prevalence",
"Prevalence_num" ))

write_delim(df_combined, snakemake@output[[1]], delim="\t")
