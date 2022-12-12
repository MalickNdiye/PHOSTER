library(tidyverse)

get_mtdata<- function(mtdata_path="./Genomes_database_metadata_MN.txt", new, corr="Locus_tag", df_corr){
  mtdata<-read.csv(mtdata_path, header = T, sep="\t")
  lt<- pull(mtdata, new)
  names(lt)<- pull(mtdata, corr)
  nt<- lt[df_corr]
  
  return(nt)
}

movMeanCirc <- function(depths, window = 100, focus = 1){
  ####
  #DESCRIPTION
  # function to compute sliding window average depth for circular chromosome
  # around a given focus base pair
  ##
  #ARGUMENTS
  # depths: vector of integers representing read depth
  # window: nr of bp before and after focus to include in average. Defaults to 500
  # focus:  index integer indicating around which bp to compute average. Defaults to 1
  ##
  #SET-UP
  # 1. define linear before-after index around focus with given window size
  # 2. find real, circular before and after index
  # 3. if before <= end: mean of values within window; else outside window
  ####
  
  # max linear index value
  linear_end <- length(depths)
  # 1. direct, linear index values
  index_left <- focus - window
  index_right <- focus + window
  # 2. real, circular index values
  index_before <- ifelse(index_left >= 1,
                         yes = index_left,
                         no = ifelse((linear_end + index_left)>=1, (linear_end + index_left),
                                     linear_end-1))
  index_after <- ifelse(index_right <= linear_end,
                        yes = index_right,
                        no = index_right - linear_end)
  # 3. mean sliding window
  res <- ifelse(index_before <= index_after,
                yes = mean(depths[index_before:index_after]),
                no = mean(depths[-((index_after + 1):(index_before - 1))]))
  # return result
  return(res)
}

get_gen_col<- function(genera, rel_ab, suppl=vector(), threshold=1, add_pal="Greys", unclass= "unclassified"){
  cols<- brewer.pal(11, "Spectral")
  
  minus<- paste("other<", threshold,"%", sep="")
  
  un<- c(minus,
         unclass)
  un_col<- c("grey2",
             "azure3")
  names(un_col)<- c(minus, unclass)
  
  levs= c("Lactobacillus",
          "Firm4",
          "Firm5",
          "Bifidobacterium",
          "Commensalibacter",
          "Bombella",
          "Acetobacteraceae",
          "Snodgrassella",
          "Gilliamella",
          "Bartonella",
          "Frischella",
          "Apis",
          "GB_VCs")
  
  col_list<- c(
    "Lactobacillus"=cols[10],
    "Firm4"="mediumpurple",
    "Firm5"=cols[10],
    "Bifidobacterium"=cols[8],
    "Commensalibacter"="tomato3",
    "Bombella"="brown",
    "Acetobacteraceae"="tomato3",
    "Snodgrassella"=cols[6],
    "Gilliamella" =cols[7],
    "Bartonella"=cols[3],
    "Frischella"=cols[5],
    "Apis"=cols[1],
    "GB_VCs"="turquoise"
  )
  
  other<- rel_ab<threshold
  genera[other]<-minus
  
  unrep_gen<- genera[!(genera %in% levs)]
  names(unrep_gen)<- unrep_gen
  unrep_gen<- unrep_gen[!(unrep_gen==minus)]
  unrep_gen<- unrep_gen[!(unrep_gen==unclass)]
  un_n<- length(unique(unrep_gen))
  
  
  new_col<-  colorRampPalette(brewer.pal(9, add_pal))(un_n)
  names(new_col)=unique(unrep_gen)
  
  col_list<- append(new_col, col_list)
  col_list<- append(un_col, col_list)
  col_list<- col_list[names(col_list) %in% genera]
  show_col(col_list)
  
  new_levs<- append(unique(unrep_gen),levs)
  new_levs<- append(un,new_levs)
  new_levs<- new_levs[new_levs %in% genera]
  
  return(list(new_levs, col_list))
}

save_plts<- function(plt, sample, out_dir){
  strain<-unique(plt$strain)
  print(unique(plt$strain))
  
  cov_hist<- ggplot(plt, aes(x=cov))+
    geom_density()+
    ggtitle(paste(unique(strain), unique(plt$SDP), sep = "\n")) +
    geom_vline(data=plt, aes(xintercept=q10), col="blue")+
    geom_vline(data=plt, aes(xintercept=q90), col="blue")+
    geom_vline(data=plt, aes(xintercept=avg), col="green")+
    theme_classic()
  
  file<- paste(out_dir, strain, "/", strain, "_", sample, "_covHist.pdf", sep="")
  dir<- paste(out_dir, strain, "/", sep="")
  
  if (file.exists(dir)){} 
  else {
    dir.create(file.path(dir))
  }
  
  pdf(file, height = 8.27, width=11.69)
  print(cov_hist)
  dev.off()
}
