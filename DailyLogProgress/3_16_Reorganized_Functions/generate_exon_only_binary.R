generate_exon_only_binary <- function(binary_folder_path, result_folder_name, custom_gene = NULL, file_pattern = "_binary.csv"){
  print(paste(Sys.time(), ": processing files.")); flush.console()
  init_files = list.files(binary_folder_path, pattern = file_pattern)
  dir.create(file.path(binary_folder_path, result_folder_name), showWarnings = FALSE)
  files <- init_files
  skipped_gene = NULL
  if(length(custom_gene) > 0){
    files = NULL
    for(i in 1:length(files)){
      temp_gene_name <- strsplit(init_files[i], split = file_pattern)[[1]]
      if(temp_gene_name %in% custom_gene){processed_file <- c(processed_file, init_files[i])}
    }
  }

  for(i in 1:length(files)){
  #for(i in 53:53){
    gene = strsplit(files[i], split = file_pattern)[[1]]
    print(paste("Processing Files for Gene #", i,":", gene))
    file <- files[i]
    df <- read.csv(file.path(binary_folder_path, file))
    isoform <- names(df)[which(df[1, ] == "R")]
    remain_binary <- df[, which(df[1, ] == "R")]

    if(length(dim(remain_binary)) != 0){
      if(dim(remain_binary)[2] != 0){
        remain_binary <- remain_binary[-1, ]
        remain_binary_mat <- as.matrix(remain_binary)
        remain_binary_mat <- apply(remain_binary_mat, 2, as.numeric)
        result <- isoform_group(remain_binary_mat, isoform)
        file_name <- paste(gene, "_exon_only_binary.csv", sep = "")
        d <- data.frame(result$exon_binary)
        group_index <- result$index
        names(d) <- isoform
        d <- rbind(group_index, d)
        write.csv(d, file.path(binary_folder_path, result_folder_name, file_name))
      }
      else{
        skipped_gene <- c(skipped_gene, gene)
      }
    }

    else{
      skipped_gene <- c(skipped_gene, gene)
    }
  }

  write(skipped_gene, file.path(binary_folder_path, result_folder_name, "Skipped_Genes.txt"))
  print(paste(Sys.time(), ": done.")); flush.console()
  print(paste("Skipped Genes:", skipped_gene))

}
