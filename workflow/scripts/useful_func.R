library(scales)

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

get_gen_col<- function(genera, rel_ab=vector(), threshold=1, add_pal="Greys", unclass= "unclassified"){
  cols<- brewer.pal(11, "Spectral")

  minus<- paste("other<", threshold,"%", sep="")
  un<- c(minus,
         unclass)
  un_col<- c("grey2",
             "azure3")
  names(un_col)<- c(minus, unclass)

  levs= c("Lactobacillus",
          "Firm4",
          "Bombilactobacillus",
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
    "Bombilactobacillus"="mediumpurple",
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

  unrep_gen<- genera[!(genera %in% levs)]

  if(length(rel_ab) >0){
    other<- rel_ab<threshold
    genera[other]<-minus
    names(unrep_gen)<- unrep_gen
    unrep_gen<- unrep_gen[!(unrep_gen==minus)]
    unrep_gen<- unrep_gen[!(unrep_gen==unclass)]
    un_n<- length(unique(unrep_gen))
  } else{
    un_n<- length(unique(unrep_gen))
  }

  new_col<-  colorRampPalette(brewer.pal(9, add_pal))(un_n)
  names(new_col)=unique(unrep_gen)

  col_list<- append(new_col, col_list)
  col_list<- append(un_col, col_list)
  col_list<- col_list[names(col_list) %in% genera]

  new_levs<- append(unique(unrep_gen),levs)
  new_levs<- append(un,new_levs)
  new_levs<- new_levs[new_levs %in% genera]

  return(list(new_levs, col_list))
}

get_gen_col_abs<- function(genera, add_pal="Greys", unclass= "unclassified"){
  cols<- brewer.pal(11, "Spectral")

  levs= c("Lactobacillus",
          "Firm4",
          "Bombilactobacillus",
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
    "Bombilactobacillus"="mediumpurple",
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

  un_col<- c("azure3")
  names(un_col)<-unclass

  unrep_gen<- genera[!(genera %in% levs)]
  un_n<- length(unique(unrep_gen))

  new_col<-  colorRampPalette(brewer.pal(9, add_pal))(un_n)
  names(new_col)=unique(unrep_gen)

  col_list<- append(new_col, col_list)
  col_list<- append(un_col, col_list)
  col_list<- col_list[names(col_list) %in% genera]

  new_levs<- append(unique(unrep_gen),levs)
  new_levs<- append(unclass,new_levs)
  new_levs<- new_levs[new_levs %in% genera]

  return(list(new_levs, col_list))
}


get_core<- function(){
  core_list<- c("Lactobacillus",
                "Firm4",
                "Bombilactobacillus",
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
  return(core_list)
}

get_core_cols<- function(){
  cols<- brewer.pal(11, "Spectral")
  col_list<- c(
    "Lactobacillus"=cols[10],
    "Firm4"="mediumpurple",
    "Bombilactobacillus"="mediumpurple",
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

  return(col_list)
}

get_species_cols<- function(species){

  genus_pal<- c("Lactobacillus"="Blues",
                "Firm4"="Purples",
                "Bombilactobacillus"="Purples",
                "Firm5"="Blues",
                "Bifidobacterium"="Set3",
                "Snodgrassella"="Oranges",
                "Gilliamella" ="PiYG",
                "Bartonella"="PuRd",
                "Friscella"="wheat",
                "Bombella"="OrRd",
                "Commensalibacter"="tomato")

  sp_col_list=c()
  for( i in names(genus_pal)){
    sp_l=unique(species[grepl(paste("\\b", i, "\\b", sep=""), species)])
    print(sp_l)
    n= length(sp_l)
    if(n >1 ){
      pal=brewer.pal(n, genus_pal[i])
      names(pal)=sp_l
      sp_col_list<- c(sp_col_list, pal)
    } else{
      pal= genus_pal[i]
      names(pal)=sp_l
      sp_col_list<- c(sp_col_list, pal)
    }
  }

  return(sp_col_list)
}
