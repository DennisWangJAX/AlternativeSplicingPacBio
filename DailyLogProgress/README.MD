Log Folder


2/19:
Added debug for exceptions:
1. Correct file input @Line 49
2. Start and end judgement in the creation of mat_gencode variable @Line 154, 155
3. Deletion of an argument (last one in vector false) in pre-filtering @Line 210

2/23:
1. Added the visualize_cluster function
2. Modified the intron_filter_by_GRanges @line 156 - 158 for the exception of having an output of vector (out_disjoin_filtered)
