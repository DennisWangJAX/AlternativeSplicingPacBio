Plotting .R

quick plot:
gff_name = "RIPK1.gtf"
results <- intron_filter_by_GRanges(import.gff(gff_name))
group_result <- isoform_group(results$remain_binary, results$remain)
#visualize_cluster(gff_name, group_result$group_index, results$remain)
visualize(gff_name, 1, transcript_name = group_result$isoform_group[[1]])
