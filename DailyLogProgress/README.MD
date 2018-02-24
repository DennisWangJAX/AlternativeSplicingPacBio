Log Folder


2/19:
Added debug for exceptions:
1. Correct file input @Line 49
2. Start and end judgement in the creation of mat_gencode variable @Line 154, 155
3. Deletion of an argument (last one in vector false) in pre-filtering @Line 210

2/23:
1. Added the visualize_cluster function
2. Modified the intron_filter_by_GRanges @line 156 - 158 for the exception of having an output of vector (out_disjoin_filtered)
3. Fixed bug in intron_filter_by_GRanges @line 36 - 38. Turn the chr value into numeric to fix the [coverage] function. A bug in IsoView as well. (Arbitrary Selection) Original Error: invalid class “GRanges” object: 'mcols(x)' is not parallel to 'x'

2/24: 
All for updates in intron_filter_by_GRanges function
1. Fixed exception at line 36 - 38. Provide the exception for X chromosome as for 23 numeric value.
2. Add if command starting at line 167. Fixed the exception when ncol(out_disjoin_filtered) == 0.
3. Change in the output, for everything that contains !test when test == NULL
