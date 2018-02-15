intron_filter <- function(gff){
  source("https://bioconductor.org/biocLite.R")
  require("stringr"); require("rtracklayer"); require("GenomicRanges")

  search_pat <- function(pattern, x){
    result <- NULL
    for (i in 1:(length(x)-length(pattern))){
      compare <- x[i:(i+length(pattern)-1)]
      if (identical(pattern,compare)){
        result <- c(result,i)
      }
    }

    result
  }

  ## Data Cleaning function
  find_retention <- function(input_binary_mat){
    test_output <- logical(ncol(input_binary_mat)) == F
    num_zero <- 1
    jud <- TRUE

    while (jud & num_zero < (nrow(input_binary_mat)-1)){
      pat <- c(1, rep(0, num_zero), 1)
      pat_1 <- rep(1, length(pat))
      for (i in 1:ncol(input_binary_mat)){
        temp <- input_binary_mat[, i]
        index_list <- search_pat(pat, temp)
        for (index in index_list){
          for (j in 1:ncol(input_binary_mat)){
            if(j != i){
              if(identical(input_binary_mat[, j][index:(index+length(pat)-1)],pat_1)){
                test_output[j] <- FALSE
              }
            }
          }
        }
      }
      num_zero <- num_zero+1
    }
    test_output
  }

  ##Prepare gtf file from hg38
  print(paste(Sys.time(), ": prepare gtf file from hg38")); flush.console()
  if(!file.exists("Gencode_proteinCoding.gtf")){
    library(AnnotationHub)
    ah = AnnotationHub()
    ah = subset(ah,species=="Homo sapiens")
    qhs <- query(ah, c("Ensembl", "gene", "annotation", "grch38"))
    gtf <- qhs[["AH51014"]] #Homo_sapiens.GRCh38.85.gtf
    export(gtf[ gtf$transcript_biotype %in% "protein_coding", ], "Gencode_proteinCoding.gtf")
  }

  gencode_gtf <- import.gff("Gencode_proteinCoding.gtf")
  gr = import.gff(gff)
  gencode_gtf_gene <- gencode_gtf[gencode_gtf$gene_name %in% gr$gffcompare_gene_name]
  gencode_exon = gencode_gtf_gene[gencode_gtf_gene$type == "exon"]
  gencode_base = reduce(gencode_exon)
  split_name <- unlist(strsplit(gff, split = "\\."))
  export_gtf_name <- paste("gencode_", split_name[1], ".gtf", sep = "")
  export(gencode_exon, export_gtf_name)
  print(paste(Sys.time(), "Compare gtf file from Gencode AH51014 is saved to:", export_gtf_name)); flush.console()

  gr_exon = gr[gr$type=="exon"]
  gr_list_input = split(gr_exon, gr_exon$transcript_id)
  isoform_names <- names(gr_list_input)
  isoform_count = length(unique(gr_exon$transcript_id))





  populateMatrix = function(query, subject){
    info = findOverlaps(query, subject)
    ind = 1:length(subject)
    return((ind %in% subjectHits(info)))
  }

  print(paste(Sys.time(), ": preprocessing")); flush.console()
  cov = base::as.vector(GenomicRanges::coverage(gr_exon)[[1]])
  extract = which(cov>0)
  extract_high = which(cov>isoform_count*0.1)
  chr = unique(seqnames(gr_list_input[[1]]))
  gr_extract = reduce(GRanges(chr, range=IRanges(start=extract, end=extract)))
  gr_extract_high = reduce(GRanges(seqnames=chr, range=IRanges(extract_high, extract_high)))

  gr_base = reduce(gr_exon)
  gr_subject = gr_extract_high
  gr_intron = GenomicRanges::setdiff(GRanges(seqname=chr, ranges = IRanges(start(gr_base), end(gr_base))), gr_subject)
  gr_exons = gr_subject; gr_exons$type = "exon"
  gr_introns = gr_intron;
  if (length(gr_introns) > 0) { gr_introns$type="intron" }
  gr_tract = c(gr_exons, gr_introns)
  gr_tract = gr_tract[order(start(gr_tract))]


  new1 <- GRanges(seqnames = chr, range = IRanges(start(gencode_exon), end(gencode_exon)))
  new2 <- GRanges(seqnames = chr, range = IRanges(start(gr_exon), end(gr_exon)))
  recombined_gencode_gr <- c(new1, new2)
  recombined_gencode_gr <- recombined_gencode_gr[order(start(recombined_gencode_gr))]
  recombined_gr_base <- reduce(recombined_gencode_gr)

  reduce_intron_start = NULL
  reduce_intron_end = NULL
  for (i in 1:(length(recombined_gr_base)-1)){
    reduce_intron_start <- c(reduce_intron_start, end(recombined_gr_base)[i]+1)
    reduce_intron_end <- c(reduce_intron_end, start(recombined_gr_base)[i+1]-1)
  }
  reduce_intron = GRanges(seqname=chr, ranges = IRanges(reduce_intron_start, reduce_intron_end))
  gr_combined = c(recombined_gr_base, reduce_intron)
  gr_combined = gr_combined[order(start(gr_combined))]


  piece = c(disjoin(recombined_gencode_gr), reduce_intron)
  piece = piece[order(start(piece))]
  piece_width = width(piece)
  piece_start = start(piece)
  piece_end = end(piece)


  print(paste(Sys.time(), ": computing similarities of exons and introns")); flush.console()
  stick = NULL
  for(i in 1:length(gr_combined)){stick=c(stick, seq(start(gr_combined[i]), end(gr_combined[i])))} #changed
  stick = GRanges(seqnames=chr, ranges =IRanges(start=stick, end=stick))
  mat = matrix(unlist(lapply(gr_list_input, populateMatrix, stick)), ncol=length(stick), byrow=T)
  mat = t(mat)
  mat_gencode = rep(FALSE, nrow(mat))

  for (i in (1:length(gencode_base))){
    gencode_start <- start(gencode_base)[i] - start(gr_base)[1]
    gencode_end <- gencode_start + width(gencode_base)[i]
    print(c(gencode_start, gencode_end))
    mat_gencode[gencode_start:gencode_end] = T
  }
  out_disjoin <- matrix(NA, ncol = ncol(mat), nrow = length(piece))
  gencode_out_disjoin <- matrix(NA, ncol= 1, length(piece))
  colnames(out) <- colnames(mat)
  rownames(out) <- colnames(mat)


  # Based on Disjoin
  for (i in 1:(ncol(out_disjoin))){
    a <- mat[, i]
    for(j in 1:(nrow(out_disjoin))){
      tract_e <- sum(width(piece)[1:j])
      tract_s <- tract_e - width(piece)[j] + 1
      compare_region <- a[tract_s:tract_e]
      val <- sum(compare_region, na.rm = TRUE)/(width(piece)[j])
      out_disjoin[j, i] <- val
    }
  }
  # Separation Line


  # Calculating for gencode standard set
  for(j in 1:(nrow(gencode_out_disjoin))){
    tract_e <- sum(width(piece)[1:j])
    tract_s <- tract_e - width(piece)[j] + 1
    compare_region <- mat_gencode[tract_s:tract_e]
    val <- sum(compare_region, na.rm = TRUE)/(width(piece)[j])
    if(val<1){val = 0}
    gencode_out_disjoin[j, 1] <- val
  }
  # Separation Line


  # Getting rid of isoforms with potential new exons
  gencode_test <- find_retention(cbind(gencode_out_disjoin, out_disjoin[,1]))

  # Separation Line


  test <- find_retention(out_disjoin)

  for (i in 1:length(isoform_names)){
    if (gr[gr$transcript_id == isoform_names[i]]$matchAnnot_gene[1] != gr[gr$transcript_id == isoform_names[i]]$gffcompare_gene_name[1]){
      test[i] <- FALSE
    }

  }


  print(paste(Sys.time(), ": exporting two files...")); flush.console()
  export(gr[gr$transcript_id %in% isoform_names[test]], paste(split_name[1],"_remain2.0",".",split_name[2], sep=""))
  export(gr[gr$transcript_id %in% isoform_names[test==F]], paste(split_name[1],"_filtered2.0",".",split_name[2], sep=""))
  list("remain" = isoform_names[test], "filtered" = isoform_names[test==F])
}
