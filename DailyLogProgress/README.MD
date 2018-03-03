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

2/26:
1. Added the outputs (csv) for the binary forms of isoforms. Each is characterized by a leading sign. "P" for prefiltered, "F" for filtered, "R" for remain.

2/27:
1. Added variable (repetition_test) in intron_filter_by_gff to save running time if the files already exist.

2/28:
1. Added new function - generate_exon_only_binary.R
2. Sample: generate_exon_only_binary(binary_folder_path = "/Users/c-wangq/analysis/SplicingPlots/binary_results/",
                          result_folder_name = "grouped_remain_binary_results",
                          custom_gene = NULL)
3. Transform all the preliminary binary forms into exon only binary forms through self-learning. 
4. First row contains unique group number with Manhattan Distance
5. Column indicates isoform names.
6. Added new function - heatmap_exon_binary.R
  plot with the input file name "*_exon_only_binary.csv"
7. Add Line 24 - 27 in isoform_group.R to fix exception when the last row is an exon.

3/2:
1. Add autogenerate_plots.R to generate new plots in a directory
2. Modified the pheatmap functions. heatmap_coverage.R
