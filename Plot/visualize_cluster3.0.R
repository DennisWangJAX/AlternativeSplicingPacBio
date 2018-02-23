visualize_cluster = function(gff, clust_count, transcript_name, out_png = "Visualize.png", customGrouping=F, grouping=NULL, random_seed=1337){
  require("stringr"); require("rtracklayer"); require("GenomicRanges")
  gr = import.gff(gff)
  gr_exon = gr[gr$type=="exon"]
  gr_exon = gr_exon[which(gr_exon$transcript_id %in% transcript_name)]
  gr_list_input = split(gr_exon, gr_exon$transcript_id)
  isoform_count = length(unique(gr_exon$transcript_id))
  
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
  num_bins = length(gr_intron) + length(gr_subject)
  gr_exons = gr_subject; gr_exons$type = "exon"
  gr_introns = gr_intron; 
  if (length(gr_introns) > 0) { gr_introns$type="intron" }
  gr_tract = c(gr_exons, gr_introns)
  gr_tract = gr_tract[order(start(gr_tract))]
  
  if(customGrouping==F){
    print(paste(Sys.time(), ": computing distance")); flush.console()
    stick = NULL
    for(i in 1:length(gr_base)){stick=c(stick, seq(start(gr_base[i]), end(gr_base[i])))}
    stick = GRanges(seqnames=chr, ranges =IRanges(start=stick, end=stick))
    mat = matrix(unlist(lapply(gr_list_input, populateMatrix, stick)), ncol=length(stick), byrow=T)
    mat = t(mat)
    out <- matrix(NA, ncol = ncol(mat), nrow = ncol(mat))
    colnames(out) <- colnames(mat)
    rownames(out) <- colnames(mat)
    for (i in 1:(ncol(mat))) {
      for (j in 1:(ncol(mat))) {
        a <- mat[, i]
        b <- mat[, j]
        out[i, j] <- (sum(b, na.rm = TRUE) + sum(a, na.rm = TRUE) - 
                        2 * sum(a & b, na.rm = TRUE))/(sum(a, na.rm = TRUE) + 
                                                         sum(b, na.rm = TRUE) - sum(a & b, na.rm = TRUE))
      }
    }
    distance=out
    set.seed(random_seed)
    clust = kmeans(distance, clust_count)
    cluster = clust$cluster
    indx = order(cluster)
  }
  else{
    cluster=grouping
    indx = order(cluster)
  }
  
  print(paste(Sys.time(), ": plotting")); flush.console()
  png(out_png, width = max(800,num_bins*80), height=max(length(indx)*15,600))
  par(mar=rep(5,4))
  plot(c(0, num_bins), c(0, length(indx)), ty="n", xaxt="n", yaxt="n", ylab="", xlab="", main="Transcript Isoform Visualization")
  rect(xleft=0:(num_bins-1), xright=1:num_bins, ytop=0, ybot=-1) #, col=(gr_tract$type=="exon")+2)
  text(x=((0:(length(gr_tract)-1)+1:length(gr_tract))/2), y = -0.1, labels = width(gr_tract))
  count=1
  widths = NULL
  for(j in indx){
    gr_subject= gr_list_input[[j]]
    ol = findOverlaps(gr_subject, gr_tract)
    qh = queryHits(ol); sh = subjectHits(ol)
    lefts = by(sh, qh, min)
    rights = by(sh, qh, max)
    ids = names(lefts)
    start = (start(gr_subject)[as.numeric(ids)] - start(gr_tract)[lefts])/width(gr_tract)[lefts]+lefts-1
    end = rights - ( (end(gr_tract)[rights]) - end(gr_subject)[as.numeric(ids)])/width(gr_tract)[rights]
    rect(xleft = start, xright=end, ytop = count, ybot = count-1, col = "blue")
    count = count+1
    widths = c(widths, sum(width(gr_subject)))
  }
  par(xpd=NA)
  nn <- names(gr_list_input)
  text(x = -(num_bins)/15, y = 0:(length(indx))+0.5, pos=3, cex=0.8, labels= c(strtrim(nn[indx], 15), "Name"))
  text(x = num_bins*16/15, y = 0:(length(indx)+1)-0.2, pos=4, cex=0.8, labels=c(sum(width(gr_tract)), widths, "Length"))
  dev.off()
  print(paste(Sys.time(), ": done")); flush.console()
  
  for (i in 1:clust_count) {
    print(paste("group", i))
    print(nn[which(cluster==i)])
  }
}