sort_nest<- function(matrix, rows=T, cols=T){
  rowO<- ifelse(rows,list(order(rowSums(matrix), decreasing = T)),list(rownames(matrix)))
  rowO<- unlist(rowO)
  
  
  colO<- ifelse(cols,
                list(order(colSums(matrix), decreasing = T)),
                list(colnames(matrix)))
  colO<-unlist(colO)
  
  sorted_mat<- matrix[rowO, colO]
  
  return(sorted_mat)
}

create_null_random<- function(mat){
  # given a binary matrix, this function returns a matrix with the same numbers of 1 but randomly distributed
  
  # get number of 1
  n_1<- sum(mat)
  
  # create null matrix
  null_mat<- matrix(0, nrow=nrow(mat), ncol=ncol(mat))
  
  # get random positions
  random_pos<- sample(1:(nrow(mat)*ncol(mat)), n_1)
  
  # fill null matrix
  null_mat[random_pos]<- 1
  
  # return null matrix
  return(null_mat)
}

create_null_prob<- function(mat){
  # given a binary matrix,return a matrix where the probability of having a 1 in position i,j is 1/2((ki/P)+(dj/H)) 
  # where ki is the number of 1 in row i, P is the total number of column in the matrix, dj is the number of 1 in column j and H is the total number of rows in the matrix
  
  # get number of 1 in each row
  n_1_row<- rowSums(mat)
  # get number of 1 in each column
  n_1_col<- colSums(mat)
  # get total number of 1
  n_1<- sum(mat)
  # get total number of rows
  n_row<- nrow(mat)
  # get total number of columns
  n_col<- ncol(mat)
  # create null matrix
  null_mat<- matrix(0, nrow=nrow(mat), ncol=ncol(mat))
  
  # fill null matrix
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      prob<- 1/2*((n_1_row[i]/n_col)+(n_1_col[j]/n_row))
      null_mat[i,j]<- rbinom(1, 1, prob)
    }
  }
  
  # return null matrix
  return(null_mat)
}



create_perfect_nest<- function(mat){
  # given a binary matrix, this function returns a perfectly nested matrix
  # it does so by "pushing" all the 1s of mat to the leftmost position that is available in the row and is not already occupied by a 1
  
  # create null matrix
  null_mat<- matrix(0, nrow=nrow(mat), ncol=ncol(mat))

  mat2<- sort_nest(mat)
  
  # fill null matrix
  for(i in 1:nrow(mat2)){
    for(j in 1:ncol(mat2)){
      if(mat2[i,j]==1){
        # find leftmost position that is available and is not already occupied by a 1
        pos<- which(null_mat[i,]==0)[1]
        
        # fill null matrix
        null_mat[i,pos]<- 1
      }
    }
  }
  
  return(null_mat)
}

measure_nestedness<- function(perfect_mat, mat){
  # given a perfectly nested matrix and a binary matrix, this function returns the nestedness of the binary matrix
  # nestedness is defined as the number of incongruent 1s and 0s in the binary matrix compared to the perfectly nested matrix
  
  # get number of incongruent 1s
  n_incongruent_1<- sum((perfect_mat==0) & (mat==1))
  # get number of incongruent 0s
  n_incongruent_0<- sum((perfect_mat==1) & (mat==0))
  
  # return nestedness
  return((n_incongruent_1+n_incongruent_0)/2)
}

calculate_p_value <- function(observed_avg, null_avg, null_sd, ci=c(-1.96, 1.96)) {
  # Calculate standard error of null normal distribution
  se <- null_sd / sqrt(length(ci))
  
  # Calculate z-score
  z_score <- (observed_avg - null_avg) / se
  
  # Calculate p-value for the upper tail
  p_upper <- 1 - pnorm(abs(z_score))
  
  # Calculate p-value for the lower tail
  p_lower <- pnorm(-abs(z_score))
  
  p_value <- p_upper + p_lower
  
  return(p_value)
}



TestNest<- function(mat, iter=10){

  print("matrix:")
  print(mat)
  
  if(length(mat)<=1){
    df<- data.frame("NODF"=0,
                    "NullRandom"=0,
                    "SDRandom"=0,
                    "NullProb"=0,
                    "SDProb" =0,
                    "p_value"=0)
    return(df)
  }
  
  random_nest<- vector(length=iter)
  prob_nest<- vector(length=iter)
  
  for(i in 1:iter){
   
    
    # create null matrices
    random_mat<- sort_nest(create_null_random(mat))
    prob_mat<- sort_nest(create_null_prob(mat))
    
    # get nestedness
    random_nest[i]<- nestednodf(random_mat)$statistic["NODF"]/100
    prob_nest[i]<-nestednodf(prob_mat)$statistic["NODF"]/100
  }
  
  avg_ranom_nestedness<- mean(random_nest)
  avg_prob_nestedness<- mean(prob_nest)
  CI_random= sd(random_nest)
  CI_prob= sd(prob_nest)
  
  original_nestedness<- nestednodf(sort_nest(mat))$statistic["NODF"]/100
  original_nestedness<- unname(original_nestedness)
  
  df<-data.frame("NODF"=original_nestedness,
                    "NullRandom"=avg_ranom_nestedness,
                    "SDRandom"=CI_random,
                    "NullProb"=avg_prob_nestedness,
                    "SDProb" =CI_prob)
  
  df<- df %>% mutate("p_value"=calculate_p_value(NODF, NullProb, SDProb, ci=c(-1.96, 1.96)))
  return(df)
}

