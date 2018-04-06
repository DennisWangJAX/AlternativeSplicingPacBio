simplify_gff <- function(gff, custom_gene = NULL, output_folder = "simplified_gtf"){
  source("https://bioconductor.org/biocLite.R")
  require("stringr"); require("rtracklayer"); require("GenomicRanges"); require("progress"); require("AnnotationHub")
  annot_gr = import.gff(gff)
  gene_names = unique(annot_gr[!is.na(annot_gr$gffcompare_gene_name)]$gffcompare_gene_name)
  count_skipped = 0
  dir.create(file.path(getwd(), "processed_gtf_files"), showWarnings = FALSE)
  output_gr_list <- NULL
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

      if (repetition_test & transcript_count > 3){
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
    for(i in 1:(length(gene_names))){
      #for(i in 1:12){
      gene = gene_names[i]
      analyzed_gene <- c(analyzed_gene, gene)
      if(gene %in% custom_gene){
        count <- count + 1
        print(paste("Analyzing Gene # ", i=count,": ", gene, sep = "" ))
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

        if (repetition_test & transcript_count > 3){
          print(paste("Match Requirements: Computing for Gene - ", gene, sep = "" ))
          temp_data <- intron_filter_by_GRanges(gr_temp)
          if(length(temp_data$remain) <= 1){
            print("Skipped Isoform Group as the number of remaining isoforms is <= 1.")
            temp_reduce <- reduce(temp_data$all_exon[temp_data$all_exon$transcript_id == temp_data$remain])
            mcols(temp_reduce) <- cbind(paste(transcript_id, "group", 1, sep = "."), "exon", transcript_id, gene)
            names(elementMetadata(temp_reduce)) <- c("transcript_id", "type", "gffcompare_gene_id", "gffcompare_gene_name")
            temp_simplified_gr <- temp_reduce
            output_file_name <- paste(gene, "_simplified.gtf", sep = "")
            export(temp_simplified_gr, file.path(getwd(), output_folder, output_file_name))
            output_gr_list <- c(output_gr_list, import.gff(output_file_name))
            }
          if(length(temp_data$remain) > 1){
            isoform_group_result <- isoform_group(temp_data$remain_binary, temp_data$remain)
            temp_simplified_gr <- simplify_GRanges(temp_data$all_exon, isoform_group_result$group, gene)
            output_gr_list <- c(output_gr_list, temp_simplified_gr)
            
          }
          }

        }

      }
    }


  #print(paste("Total Skipped Genes: ", count_skipped))
  output_gr_list <- unlist(GRangesList(output_gr_list))
  if(length(output_gr_list) > 0){export(output_gr_list, paste("simplified_", gff, sep = ""))}
  list("output" = output_gr_list, "computed_genes" = analyzed_gene)



}
