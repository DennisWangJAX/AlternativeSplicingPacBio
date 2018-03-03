heatmap_coverage <- function(gr_input, output_file_path = getwd(), isoform_group_index){
  gene_name <- unique(gr_input$gffcompare_gene_name)[1]
  dir.create(file.path(output_file_path, "output_heatmaps_2.0/"), showWarnings = F)
  gr_exon <- gr_input[gr_input$type == "exon"]
  gr_list_input = split(gr_exon, gr_exon$transcript_id)
  isoform_names <- names(gr_list_input)
  isoform_count = length(unique(gr_exon$transcript_id))

  chr = unique(seqnames(gr_list_input[[1]]))
  chr_num = as.numeric(strsplit(toString(chr), split = "chr")[[1]][2])
  if(is.na(chr_num)){chr_num = 23}
  cov = base::as.vector(GenomicRanges::coverage(gr_exon)[[1]])
  extract = which(cov>0)
  extract_high = which(cov>isoform_count*0.1)

  gr_extract = reduce(GRanges(chr, range=IRanges(start=extract, end=extract)))
  gr_extract_high = reduce(GRanges(seqnames=chr, range=IRanges(extract_high, extract_high)))

  gr_base = reduce(gr_exon)
  stick = NULL
  for(i in 1:length(gr_base)){stick=c(stick, seq(start(gr_base[i]), end(gr_base[i])))}
  stick = GRanges(seqnames=chr, ranges =IRanges(start=stick, end=stick))
  mat = matrix(unlist(lapply(gr_list_input, populateMatrix, stick)), ncol=length(stick), byrow=T)
  mat = t(mat)
  coverage_out = matrix(NA, nrow = length(gr_base), ncol = length(gr_list_input))


  for(i in 1:length(gr_base)){
    index_end <- sum(width(gr_base)[1:i])
    index_start <- index_end - width(gr_base)[i] + 1
    wid <- width(gr_base)[i]
    for(j in 1:length(gr_list_input)){
      val <- sum(mat[, j][index_start:index_end])/wid
      coverage_out[i, j] <- val
    }
  }
  coverage_out <- t(coverage_out)
  rownames(coverage_out) <- isoform_names
  colnames(coverage_out) <- width(gr_base)

  quantity <- NULL
  final_coverage_out <- matrix(NA, nrow = length(unique(iso_index)), ncol = length(gr_base))
  for(i in 1:length(unique(iso_index))){
    if(length(which(iso_index == i)) == 1){final_coverage_out[i, ] <- coverage_out[which(iso_index == i), ]}
    else{
      for(j in 1:length(gr_base)){
        match_index <- which(iso_index == i)
        val <- sum(coverage_out[match_index, ][, j])/length(match_index)
        final_coverage_out[i, j] <- val
      }
    }
    quantity <- c(quantity, length(which(iso_index == i)))
  }

  rownames(final_coverage_out) <- paste("M = ", quantity, sep = "")
  colnames(final_coverage_out) <- width(gr_base)


  out_png <- file.path(output_file_path, "output_heatmaps_2.0", paste(gene, "_coverage_heatmap.png", sep = ""))
  myColor <- colorRampPalette(c("white", "blue", "black"))(100)
  mat <- matrix(NA, nrow = nrow(final_coverage_out), ncol = ncol(final_coverage_out))
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      if(final_coverage_out[i, j] != 0){mat[i, j] <- round(final_coverage_out[i, j], 2)}
      else{mat[i, j] <- NA}
    }
  }
  cluster_r = FALSE
  if(length(unique(isoform_group_index)) > 1){cluster_r = TRUE}
  
  pheatmap(final_coverage_out, cluster_rows = cluster_r, cluster_cols = FALSE, display_numbers = TRUE, color = myColor, main = paste("Isoform Groups of", gene), filename = out_png,
           cellwidth = 48.5, cellheight = 30, number_color = "white")
  
}
