isoform_group <- function(input_binary_mat, isoform_names){
  require(mgcv)
  if (length(isoform_names) != ncol(input_binary_mat)){print("Warning: Isoform Names Input might not match the binary matrix.")}
  exon_start_index = NULL
  exon_end_index = NULL
  index = 1
  standard_exon_binary_list = rep(0, nrow(input_binary_mat))
  while (index < nrow(input_binary_mat)){
    if (sum(input_binary_mat[index, ]) > 0){
      temp = 0
      while((index+temp) < nrow(input_binary_mat) & sum(input_binary_mat[(index+temp), ]) > 0){
        temp <- temp + 1
      }
      exon_start_index <- c(exon_start_index, index)
      exon_end_index <- c(exon_end_index, index+temp-1)
      standard_exon_binary_list[index:(index+temp-1)] = 1
      index <- index + temp
    }

    else{index <- index + 1}

  }

  if(sum(input_binary_mat[nrow(input_binary_mat), ]) > 0 & sum(input_binary_mat[nrow(input_binary_mat)-1, ]) == 0){
    exon_start_index <- c(exon_start_index, nrow(input_binary_mat))
    exon_end_index <- c(exon_end_index, nrow(input_binary_mat))
  }

  exon_binary_mat <- matrix(0, nrow = length(exon_start_index), ncol = ncol(input_binary_mat))
  for (j in (1:ncol(input_binary_mat))){
    for(i in (1:nrow(exon_binary_mat))){
      if (sum(input_binary_mat[, j][exon_start_index[i]:exon_end_index[i]]) > 0){
        exon_binary_mat[i, j] <- 1
      }
    }
  }

  category <- uniquecombs(t(exon_binary_mat))
  group_index <- attr(category, "index")
  cluster = group_index
  indx <- order(cluster)
  pre_group <- list()

  for (i in 1:length(unique(group_index))){
    pos = which(group_index==i)
    pre_group[[i]] <- isoform_names[pos]
  }

  list("index"= group_index, "exon_binary" = exon_binary_mat, "group" = pre_group)

}
