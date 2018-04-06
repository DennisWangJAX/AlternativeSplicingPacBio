group_gff <- function(gff, custom_gene = NULL){
  source("https://bioconductor.org/biocLite.R")
  require("stringr"); require("rtracklayer"); require("GenomicRanges"); require("qpcR"); require("AnnotationHub"); require("openxlsx")
  annot_gr = import.gff(gff)
  gene_names = unique(annot_gr[!is.na(annot_gr$gffcompare_gene_name)]$gffcompare_gene_name)
  count_skipped = 0
  output_matrix <- list()
  analyzed_gene <- NULL

  if (length(custom_gene) == 0){
      count = 0
      for(i in 1:(length(gene_names))){
      #for(i in 1:1000){
      gene = gene_names[i]
      count <- count + 1
      print(paste("Analyzing Gene # ", count,": ", gene, sep = "" ))
      gr_temp_inter = annot_gr[annot_gr$gffcompare_gene_name %in% gene]
      temp_gene_id = unique(gr_temp_inter$gene_id)[1]
      gr_temp = annot_gr[annot_gr$gene_id == temp_gene_id]
      temp_gr_exon = gr_temp[gr_temp$type == "exon"]
      transcript_count = length(unique(gr_temp$transcript_id))
      repetition_test = TRUE
      if(file.exists(file.path(getwd(), output_folder, paste(gene, "_simplified.gtf", sep = "")))){
        repetition_test = FALSE
        print(paste("Match Requirements: Skipped Completed Computation for Gene - ", gene, sep = "" ))
        output_gr_list <- c(output_gr_list, import.gff(file.path(getwd(), output_folder, paste(gene, "_simplified.gtf", sep = ""))))
      }

      if (repetition_test & transcript_count > 3 & length(temp_gr_exon) > 2*transcript_count & nchar(gene) < 13){
        print(paste("Match Requirements: Computing for Gene - ", gene, sep = "" ))
        analyzed_gene <- c(analyzed_gene, gene)
        temp_data <- intron_filter_by_GRanges(gr_temp)
        if(length(temp_data$remain) <= 1){
          print("Skipped Isoform Group as the number of remaining isoforms is <= 1.")
          temp_reduce <- reduce(temp_data$all_exon[temp_data$all_exon$transcript_id == temp_data$remain])
          mcols(temp_reduce) <- cbind(paste(transcript_id, "group", 1, sep = "."), "exon", transcript_id, gene)
          names(elementMetadata(temp_reduce)) <- c("transcript_id", "type", "gffcompare_gene_id", "gffcompare_gene_name")
          temp_simplified_gr <- temp_reduce
          output_file_name <- paste(gene, "_simplified.gtf", sep = "")
          export(temp_simplified_gr, file.path(getwd(), output_folder, output_file_name))
          output_gr_list <- c(output_gr_list, temp_simplified_gr)
          }
        if(length(temp_data$remain) > 1){
          isoform_group_result <- isoform_group(temp_data$remain_binary, temp_data$remain)
          temp_simplified_gr <- simplify_GRanges(temp_data$all_exon, isoform_group_result$group, gene)
          output_gr_list <- c(output_gr_list, temp_simplified_gr)

        }
        }

      }


    }


  else{
    count = 0
    for(i in 1:(length(custom_gene))){
    #for(i in 3:5){
      gene = custom_gene[i]
      print(paste("Analyzing Gene # ", i,": ", gene, sep = "" ))
      gr_temp_inter = annot_gr[annot_gr$gffcompare_gene_name %in% gene]
      temp_gene_id = unique(gr_temp_inter$gene_id)[1]
      if(!is.na(temp_gene_id)){
        gr_temp = annot_gr[annot_gr$gene_id == temp_gene_id]
        transcript_count = length(unique(gr_temp$transcript_id))
        if (transcript_count > 3){
          analyzed_gene <- c(analyzed_gene, gene)
          print(paste("Match Requirements: Computing for Gene - ", gene, sep = "" ))
          temp_data <- intron_filter_by_GRanges(gr_temp)
          if(length(temp_data$remain) > 1){
            isoform_group_result <- isoform_group(temp_data$remain_binary, temp_data$remain)
            output_matrix <- qpcR:::rbind.na(output_matrix, c(gene, t(isoform_group_result$group)))
          }
          else{
            output_matrix <- qpcR:::rbind.na(output_matrix, c(gene, temp_data$remain))
          }
        }
      }

    }

  }


  output_matrix <- output_matrix[-1,]
  rownames(output_matrix) <- output_matrix[,1]
  output_matrix <- output_matrix[,-1]
  colnames(output_matrix) <- paste("Group_", 1:ncol(output_matrix), sep="")
  write.xlsx(data.frame(output_matrix), "result.csv", row.names = TRUE)
  list("output" = output_matrix, "computed_genes" = analyzed_gene)



}
