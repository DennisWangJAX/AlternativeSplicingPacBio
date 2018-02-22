## Data Cleaning function
find_retention <- function(input_binary_mat, head_tail_test = F){
  test_output <- logical(ncol(input_binary_mat)) == F
  num_zero <- 1
  jud <- TRUE

  while (jud & num_zero < (nrow(input_binary_mat)-1)){
    pat <- c(1, rep(0, num_zero), 1)
    pat_1 <- rep(1, length(pat))
    for (i in 1:ncol(input_binary_mat)){
      temp <- input_binary_mat[, i]
      index_list <- search_pat(pat, temp)
      for (index in index_list){
        for (j in 1:ncol(input_binary_mat)){
          if(j != i){
            if(identical(input_binary_mat[, j][index:(index+length(pat)-1)],pat_1)){
              test_output[j] <- FALSE
              if(head_tail_test){
                for(k in (1:ncol(input_binary_mat))){
                  if(any(is.na(input_binary_mat[, k][index:(index+length(pat)-1)]))){test_output[k] = F}
                }
              }
            }
          }
        }
      }
    }
    num_zero <- num_zero+1
  }
  test_output
}
