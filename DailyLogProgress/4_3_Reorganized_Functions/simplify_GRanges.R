simplify_GRanges <- function(input_gr_exon, group_list, input_gene_name, output_folder = "simplified_gtf"){
  addition_gtf <- GRanges(NULL)
  transcript_id <- paste(unlist(strsplit(input_gr_exon$transcript_id[1], split = "[.]"))[1], unlist(strsplit(input_gr_exon$transcript_id[1], split = "[.]"))[2], sep = ".")
  for (i in 1:length(group_list)){
    subgroup_isoforms <- group_list[[i]]
    subgroup_exon <- input_gr_exon[input_gr_exon$transcript_id %in% subgroup_isoforms]
    temp_reduce <- reduce(subgroup_exon)
    mcols(temp_reduce) <- cbind(paste(transcript_id, "group", i, sep = "."), "exon", transcript_id, input_gene_name)
    names(elementMetadata(temp_reduce)) <- c("transcript_id", "type", "gffcompare_gene_id", "gffcompare_gene_name")
    addition_gtf <- c(addition_gtf, temp_reduce)
  }

  #names(addition_gtf) = paste("Group", 1:length(addition_gtf), sep = "")
  dir.create(file.path(getwd(), output_folder), showWarnings = FALSE)
  output_file_name <- paste(input_gene_name, "_simplified.gtf", sep = "")
  export(addition_gtf, file.path(getwd(), output_folder, output_file_name))
  addition_gtf
}
