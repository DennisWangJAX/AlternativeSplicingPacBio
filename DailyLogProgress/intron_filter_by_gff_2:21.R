intron_filter_by_gff <- function(gff){
  source("https://bioconductor.org/biocLite.R")
  require("stringr"); require("rtracklayer"); require("GenomicRanges")

  annot_gr = import.gff(gff)
  gene_names = unique(annot_gr[!is.na(annot_gr$gffcompare_gene_name)]$gffcompare_gene_name)
  stat_prefiltered = NULL
  stat_remain = NULL
  stat_filtered = NULL
  perc_prefiltered = NULL
  perc_remain = NULL
  perc_filtered = NULL
  analyzed_gene = NULL
  count_skipped = 0

  #for(i in 1:(length(gene_names))){
  for(i in 1:10){
    gene = gene_names[i]
    gr_temp_inter = annot_gr[annot_gr$gffcompare_gene_name %in% gene]
    temp_gene_id = unique(gr_temp_inter$gene_id)[1]
    gr_temp = annot_gr[annot_gr$gene_id == temp_gene_id]
    temp_gr_exon = gr_temp[gr_temp$type == "exon"]
    transcript_count = length(unique(gr_temp$transcript_id))
    if (transcript_count > 1 & length(temp_gr_exon) > 2*transcript_count ){
      print(paste("Analyzing ", gene))
      analyzed_gene <- c(analyzed_gene, gene)
      temp_data <- intron_filter_by_GRanges(gr_temp)
      stat_prefiltered <- c(stat_prefiltered, temp_data$prefiltered_new_exon_count)
      stat_remain <- c(stat_remain, temp_data$remain_count)
      stat_filtered <- c(stat_filtered, temp_data$filtered_count)
      isoform_num <- temp_data$isoform_count
      perc_prefiltered <- c(perc_prefiltered, temp_data$prefiltered_new_exon_count/isoform_num)
      perc_remain <- c(perc_remain, temp_data$remain_count/isoform_num)
      perc_filtered <- c(perc_filtered, temp_data$filtered_count/isoform_num)
    }

    else{count_skipped <- count_skipped + 1}
  }

  print(paste("Total Skipped Genes: ", count_skipped))
  list("stat_prefiltered" = stat_prefiltered, "stat_remain" = stat_remain, "stat_filtered" = stat_filtered,
       "analyzed_gene" = analyzed_gene, "percent_filtered" = perc_filtered,
       "percent_prefiltered" = perc_prefiltered, "percent_remain" = perc_remain)
}
