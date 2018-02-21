visualize_cluster = function(gff, cluster, transcript_name){
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
  num_bins = length(gr_base)
  gr_exons = gr_subject; gr_exons$type = "exon"
  gr_introns = gr_intron; 
  if (length(gr_introns) > 0) { gr_introns$type="intron" }
  
  
  indx = order(cluster)
  print(paste(Sys.time(), ": plotting")); flush.console()
  #png(out_png, width = max(800,num_bins*80), height=max(length(indx)*15,600))
  par(mar=rep(5,4))
  plot(c(0, num_bins), c(0, length(indx)), ty="n", xaxt="n", yaxt="n", ylab="", xlab="", main="Transcript Isoform Visualization")
  rect(xleft=0:(num_bins-1), xright=1:num_bins, ytop=0, ybot=-1) #, col=(gr_base$type=="exon")+2)
  text(x=((0:(length(gr_base)-1)+1:length(gr_base))/2), y = -0.2, labels = width(gr_base), cex = 0.5)
  count=1
  widths = NULL
  for(j in indx){
    gr_subject= gr_list_input[[j]]
    ol = findOverlaps(gr_subject, gr_base)
    qh = queryHits(ol); sh = subjectHits(ol)
    lefts = by(sh, qh, min)
    rights = by(sh, qh, max)
    ids = names(lefts)
    start = (start(gr_subject)[as.numeric(ids)] - start(gr_base)[lefts])/width(gr_base)[lefts]+lefts-1
    end = rights - ( (end(gr_base)[rights]) - end(gr_subject)[as.numeric(ids)])/width(gr_base)[rights]
    rect(xleft = start, xright=end, ytop = count, ybot = count-1, col = sort(cluster)[count])
    count = count+1
    widths = c(widths, sum(width(gr_subject)))
  }
  par(xpd=NA)
  nn <- names(gr_list_input)
  #text(x = -(num_bins)/15, y = 0:(length(indx))-0.5, pos=3, cex=0.8, labels= c(strtrim(nn[indx], 10), "Name"))
  #text(x = num_bins*16/15, y = 0:(length(indx)+1)-0.5, pos=4, cex=0.8, labels=c(sum(width(gr_base)), widths, "Length"))
  text(x = -(num_bins)/15 - 1.7, y = 0:(length(indx)) - 0.25, pos=3, cex=0.7, labels= c(strtrim(nn[indx], 10), "Name"))
  text(x = num_bins*16/15, y = 0:(length(indx)+1)-0.4, pos=4, cex=0.7, labels=c(sum(width(gr_base)), widths, "Length"))
  #dev.off()
  print(paste(Sys.time(), ": done")); flush.console()
  
}