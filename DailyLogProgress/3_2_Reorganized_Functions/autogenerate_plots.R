file_list <- list.files(pattern = "_exon_only_binary.csv")

for (i in 1:length(file_list)){
  temp_final <- read.csv(file_list[i])
  gene = strsplit(file_list[i], split = "_exon_only_binary.csv")[[1]]
  iso_names = names(temp_final)[-1]
  iso_index <- as.numeric(temp_final[1, ][-1])
  file_path <- file.path("/Users/c-wangq/analysis/SplicingPlots", paste(gene, "_remain2.0.gtf", sep = ""))
  temp = import.gff(file_path) 
  #gr_exon <- temp[temp$type == "exon"]
  
  #if(length(reduce(gr_exon)) != nrow(temp_final) - 1){
  #  print(gene)
    #print(i)
    #print(nrow(temp_final))
  #  print(length(reduce(temp_exon)) - nrow(temp_final))
  print(paste("Analyzing Gene #", i, ":", gene))
  if(nrow(temp_final) > 2){
    heatmap_coverage(temp, getwd(), iso_index)
  }
}
