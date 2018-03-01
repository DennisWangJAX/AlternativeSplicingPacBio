intron_filter_by_gff <- function(gff, custom_gene = NULL){
  source("https://bioconductor.org/biocLite.R")
  require("stringr"); require("rtracklayer"); require("GenomicRanges"); require("progress"); require("AnnotationHub")
  annot_gr = import.gff(gff)
  gene_names = unique(annot_gr[!is.na(annot_gr$gffcompare_gene_name)]$gffcompare_gene_name)
  stat_prefiltered = NULL
  stat_remain = NULL
  stat_filtered = NULL
  stat_isoform = NULL
  perc_prefiltered = NULL
  perc_remain = NULL
  perc_filtered = NULL
  analyzed_gene = NULL
  count_skipped = 0
  #exception <- c("PSMB9")
  #custom_gene <- custom_gene[!custom_gene %in% exception]
  dir.create(file.path(getwd(), "binary_results"), showWarnings = FALSE)

  if (length(custom_gene) == 0){
    for(i in 1:(length(gene_names))){
      #for(i in 1:10){
      gene = gene_names[i]
      gr_temp_inter = annot_gr[annot_gr$gffcompare_gene_name %in% gene]
      temp_gene_id = unique(gr_temp_inter$gene_id)[1]
      gr_temp = annot_gr[annot_gr$gene_id == temp_gene_id]
      temp_gr_exon = gr_temp[gr_temp$type == "exon"]
      transcript_count = length(unique(gr_temp$transcript_id))
      repetition_test = TRUE
      if(file.exists(paste("gencode_", gene, ".gtf", sep = ""))
         & file.exists(paste(getwd(), "/binary_results/", gene, "_binary.csv", sep = ""))
         & (file.exists(paste(gene,"_remain2.0.gtf", sep=""))
            | file.exists(paste(gene,"_potential_new_exon.gtf", sep=""))
            | file.exists(paste(gene,"_filtered2.0.gtf", sep="")))){
        repetition_test = FALSE
      }

      if (repetition_test & transcript_count > 5 & length(temp_gr_exon) > 2*transcript_count ){
        analyzed_gene <- c(analyzed_gene, gene)
        temp_data <- intron_filter_by_GRanges(gr_temp)
        stat_prefiltered <- c(stat_prefiltered, temp_data$prefiltered_new_exon_count)
        stat_remain <- c(stat_remain, temp_data$remain_count)
        stat_filtered <- c(stat_filtered, temp_data$filtered_count)
        isoform_num <- temp_data$isoform_count
        temp_data$remain
        stat_isoform <- c(stat_isoform, isoform_num)
        perc_prefiltered <- c(perc_prefiltered, temp_data$prefiltered_new_exon_count/isoform_num)
        perc_remain <- c(perc_remain, temp_data$remain_count/isoform_num)
        perc_filtered <- c(perc_filtered, temp_data$filtered_count/isoform_num)
        d <- data.frame(stat_prefiltered, perc_prefiltered, stat_filtered, perc_filtered, stat_remain, perc_remain, stat_isoform, analyzed_gene)
        names(d) <- c("Count_Prefiltered", "Percent_Prefiltered", "Count_Filtered", "Percent_Filtered", "Count_Remain", "Percent_Remain", "Count_Isoforms", "Gene_Name")
        write.csv(d, "result.csv")
        category = c(rep("P", temp_data$prefiltered_new_exon_count), rep("F", temp_data$filtered_count), rep("R", temp_data$remain_count))
        binary_data = data.frame(temp_data$prefiltered_binary, temp_data$filtered_binary, temp_data$remain_binary)
        #print(c(temp_data$prefiltered_new_exon, temp_data$filtered, temp_data$remain))
        names(binary_data) <- c(temp_data$prefiltered_new_exon, temp_data$filtered, temp_data$remain)
        binary_data <- rbind(category, binary_data)
        bin_file_name <- paste(gene, "_binary.csv", sep = "")
        write.csv(binary_data, paste(getwd(), "/binary_results/", bin_file_name, sep = ""))
      }


    }
  }

  else{
    count = 0
    for(i in 1:(length(gene_names))){
      #for(i in 1:12){
      gene = gene_names[i]
      if(gene %in% custom_gene){
        count <- count + 1
        print(paste("Analyzing Gene # ", count,": ", gene, sep = "" ))
        gr_temp_inter = annot_gr[annot_gr$gffcompare_gene_name %in% gene]
        temp_gene_id = unique(gr_temp_inter$gene_id)[1]
        gr_temp = annot_gr[annot_gr$gene_id == temp_gene_id]
        temp_gr_exon = gr_temp[gr_temp$type == "exon"]
        transcript_count = length(unique(gr_temp$transcript_id))
        repetition_test = TRUE
        if(file.exists(paste("gencode_", gene, ".gtf", sep = ""))
           & file.exists(paste(getwd(), "/binary_results/", gene, "_binary.csv", sep = ""))
           & (file.exists(paste(gene,"_remain2.0.gtf", sep=""))
              | file.exists(paste(gene,"_potential_new_exon.gtf", sep=""))
              | file.exists(paste(gene,"_filtered2.0.gtf", sep="")))){
          repetition_test = FALSE
        }

        if (repetition_test & transcript_count > 1 & length(temp_gr_exon) > 2*transcript_count ){
          analyzed_gene <- c(analyzed_gene, gene)
          temp_data <- intron_filter_by_GRanges(gr_temp)
          stat_prefiltered <- c(stat_prefiltered, temp_data$prefiltered_new_exon_count)
          stat_remain <- c(stat_remain, temp_data$remain_count)
          stat_filtered <- c(stat_filtered, temp_data$filtered_count)
          isoform_num <- temp_data$isoform_count
          stat_isoform <- c(stat_isoform, isoform_num)
          perc_prefiltered <- c(perc_prefiltered, temp_data$prefiltered_new_exon_count/isoform_num)
          perc_remain <- c(perc_remain, temp_data$remain_count/isoform_num)
          perc_filtered <- c(perc_filtered, temp_data$filtered_count/isoform_num)
          d <- data.frame(stat_prefiltered, perc_prefiltered, stat_filtered, perc_filtered, stat_remain, perc_remain, stat_isoform, analyzed_gene)
          names(d) <- c("Count_Prefiltered", "Percent_Prefiltered", "Count_Filtered", "Percent_Filtered", "Count_Remain", "Percent_Remain", "Count_Isoforms", "Gene_Name")
          #print(d)
          write.csv(d, "result.csv")
          category = c(rep("P", temp_data$prefiltered_new_exon_count), rep("F", temp_data$filtered_count), rep("R", temp_data$remain_count))
          binary_data = data.frame(temp_data$prefiltered_binary, temp_data$filtered_binary, temp_data$remain_binary)
          #print(c(temp_data$prefiltered_new_exon, temp_data$filtered, temp_data$remain))
          names(binary_data) <- c(temp_data$prefiltered_new_exon, temp_data$filtered, temp_data$remain)
          binary_data <- rbind(category, binary_data)
          bin_file_name <- paste(gene, "_binary.csv", sep = "")
          write.csv(binary_data, paste(getwd(), "/binary_results/", bin_file_name, sep = ""))
        }

      }
    }
  }

  #print(paste("Total Skipped Genes: ", count_skipped))

  list("stat_prefiltered" = stat_prefiltered, "stat_remain" = stat_remain, "stat_filtered" = stat_filtered,
       "analyzed_gene" = analyzed_gene, "percent_filtered" = perc_filtered,
       "percent_prefiltered" = perc_prefiltered, "percent_remain" = perc_remain)
}
