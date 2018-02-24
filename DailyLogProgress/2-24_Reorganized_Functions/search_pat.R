search_pat <- function(pattern, x){
  result <- NULL
  for (i in 1:(length(x)-length(pattern))){
    compare <- x[i:(i+length(pattern)-1)]
    if (identical(pattern,compare)){
      result <- c(result,i)
    }
  }

  result
}
