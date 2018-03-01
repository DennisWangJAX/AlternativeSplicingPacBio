heatmap_coverage <- function(gr_input, output_file_path = getwd()){
  gene_name <- unique(gr_input$gffcompare_gene_name)[1]
  dir.create(file.path(output_file_path, "output_heatmaps"), showWarnings = FALSE)
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

  out_png <- file.path(output_file_path, "output_heatmaps", paste(gene, "_coverage_heatmap.png", sep = ""))
  png(out_png, width = max(800,length(gr_base)*80), height=max(length(isoform_names)*15,600))
  pheatmap(coverage_out, cluster_cols = FALSE, display_numbers = TRUE)
  dev.off()
}
