quantify <- function(group_csv, quant_file){
  require("openxlsx")
  quant_data <- read.table(quant_file)
  group_info <- read.xlsx(group_csv)
  output_data <- list()
  for(i in 1:length(group_info[, 1])){
    gene = group_info[i, 1]
    print(paste("Analyzing Gene #", i, "-", gene))
    temp_no_na <- group_info[i, ][!is.na(group_info[i,])][-1]
    for(j in 1:length(temp_no_na)){
      temp_col_name <- paste(gene, j, sep = "_")
      temp_isoform_names <- unlist(strsplit(temp_no_na[j], split = ", "))
      temp_sum <- rep(0, ncol(quant_data)-1)
      for(k in 1:length(temp_isoform_names)){
        temp_transcript_id = temp_isoform_names[k]
        temp_index <- which(quant_data[,1] == temp_transcript_id)
        temp_row_data<- quant_data[temp_index,][-1]
        temp_processed_row_data <- as.numeric(as.vector(as.matrix(temp_row_data)))
        temp_sum <- rowSums(cbind(temp_sum, temp_processed_row_data), na.rm = TRUE)
      }
      output_data <- rbind(output_data, c(paste(gene, "group", j, sep = "_"), temp_sum))

    }
  }
  colnames(output_data) <- c("gene_name", as.vector(as.matrix(quant_data[1,]))[-1])
  write.csv(output_data, "quantity.csv", row.names = FALSE)
}
