heatmap_exon_binary <- function(input_file){
  c <- read.csv(input_file)
  isoform <- names(c)
  c <- c[-1, ]
  c_mat <- apply(c, 2, as.numeric)
  c_mat <- t(c_mat[, -1])
  colnames(c_mat) <- paste("E", 1:ncol(c_mat), sep = "")
  pheatmap(c_mat, cluster_cols = FALSE)
}
