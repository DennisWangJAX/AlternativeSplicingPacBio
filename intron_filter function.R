intron_filter <- function(gff){
  require("stringr"); require("rtracklayer"); require("GenomicRanges")
  gr = import.gff(gff)
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



  reduce_intron_start = NULL
  reduce_intron_end = NULL
  for (i in 1:(length(gr_base)-1)){
    reduce_intron_start <- c(reduce_intron_start, end(gr_base)[i]+1)
    reduce_intron_end <- c(reduce_intron_end, start(gr_base)[i+1]-1)
  }
  reduce_intron = GRanges(seqname=chr, ranges = IRanges(reduce_intron_start, reduce_intron_end))
  gr_combined = c(gr_base, reduce_intron)
  gr_combined = gr_combined[order(start(gr_combined))]

  piece = c(disjoin(gr_exon), reduce_intron)
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
  out <- matrix(NA, ncol = ncol(mat), nrow = length(gr_tract))
  out_binary <- matrix(NA, ncol = ncol(mat), nrow = length(gr_tract))
  out_disjoin <- matrix(NA, ncol = ncol(mat), nrow = length(piece))
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

  ## Data Cleaning
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


  test <- logical(ncol(out_disjoin)) == F
  num_zero <- 1
  jud <- TRUE

  while (jud & num_zero < (nrow(out_disjoin)-1)){
    pat <- c(1, rep(0, num_zero), 1)
    pat_1 <- rep(1, length(pat))
    initial_length <- sum(test)
    for (i in 1:ncol(out_disjoin)){
      temp <- out_disjoin[, i]
      index_list <- search_pat(pat, temp)
      for (index in index_list){
        for (j in 1:ncol(out_disjoin)){
          if(j != i){
            if(identical(out_disjoin[, j][index:(index+length(pat)-1)],pat_1)){
              test[j] <- FALSE
            }
          }
        }
      }
    }
    num_zero <- num_zero+1

  }

  for (i in 1:length(isoform_names)){
    if (gr[gr$transcript_id == isoform_names[i]]$matchAnnot_gene[1] != gr[gr$transcript_id == isoform_names[i]]$gffcompare_gene_name[1]){
      test[i] <- FALSE
    }

  }

  
  print(paste(Sys.time(), ": exporting two files...")); flush.console()
  split_name <- unlist(strsplit(gff, split = "\\."))
  export(gr[gr$transcript_id %in% isoform_names[test]], paste(split_name[1],"_remain",".",split_name[2], sep=""))
  export(gr[gr$transcript_id %in% isoform_names[test==F]], paste(split_name[1],"_filtered",".",split_name[2], sep=""))
  list("remain" = isoform_names[test], "filtered" = isoform_names[test==F])
}
