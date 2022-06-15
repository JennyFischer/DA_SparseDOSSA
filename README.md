# DA_SparseDOSSA
Differential analysis of SparseDOSSA output with different methods and first doing normalisation



I divided the data in the different analysis approaches. The resulting files which I used for my analysis are included in the files and scripts. Some findings are printed in the terminal, like for the file Netcomi_toos_whole_analysis.R some other results are save as files, like in the script NetCoMi_tools.R or Run_Rebacca.R. Due to further analysis of the outputs of for example NetCoMi_tools.R or Run_Rebacca.R it made sense to save them.
For further information about the analysis and results I included my Bachelor Thesis.
Differential abundance analysis:
separate_output.py:
separates the output of SPARSEDOSSA2. SPARSEDOSSA2 includes three tables in its output, only two are needed. Therefore the output data is splitted into the three tables. The spike-in and absolute abundance tables are saved for further analysis.
aldex2.R
has the three output files of separate_output.py as input
Output: A dataframe which shows the intersection of found significant OTUs. If empty then there are no significant OTUs which are found in absolute as well as relative abundance. At the moment it runs for the first dataset, the paths in line 11, 15 and 22 need to be exchanged
norm_wilcoxon.R
performs Wilcoxon test for dataset(at the moment dataset 1).
Output: for each input a table which include p-values calculated for each OTU by wilcoxon rank sum test.
“gene_name” has to be added to the data frame, before performing wilcoxon. There is a notation each time in the code.
The paths also need to be changed depending on where the input data is.
norm_ttest.R
TTest for dataset(at the moment dataset 1).
Output: for each input a table which include p-values calculated for each OTU by t-test. “gene_name” has to be added to the data frame, before performing wilcoxon. There is a notation each time in the code.
The paths also need to be changed depending on where the input data is.
OTU_distribution_plot.R
was performed to see the different abundances of the absolute and relative data. In barplots. Uses raw counts of abundances as input.
OTU_intersection.R
calculates the intersection of equal significant OTUs found by all tools.
 upset_plot.R
used to plot the intersections of found significant OTUs, plot 3.2 in thesis.
scatter_plot.R
was used for the plot 3.4 in the thesis.
Microbial co-occurence network:
NetCoMI_tools.R
This script is used to perform SparCC, CCLasso, SPIEC-EASI, Pearson and Spearman with NetCoMi. The output is further used for network analysis in the Netcomi_tools_whole_analysis.R
The input is the raw count matrices for absolute and relative abundance. The path need to be changed in case of using different datasets.
The output is saved as Data files for further usage.
Netcomi_tools_whole_analysis.R
This script performs the analysis of all used tools besides REBACCA. Included are the metrics and network topology analysis.
If usage for another dataset is of interest, then the path and data name should be changed in line 15(absolute abundance) and 16(relative abundance).
All the output is printed in the terminal.
Run_Rebacca.R
This script needs the additional script REBACCA_main.R in order to run.
The script performs using also the raw count matrices of absolute and relative abundance. The path also need to be changed, when using different data.
Rebacca_whole_analysis.R
Performs network analysis for REBACCA equal to the NetCoMi analysis in Netcomi_tools_whole_analysis.R
The input is the correlation estimation matrix calculated from Run_Rebacca.R
test_clr.R
This script was used to test whether there is a performance difference between absolute and relative abundance, when you use clr-transformation for relative abundance. Therefore Spearman and Pearson were performed as unnormalised absolute abundance, unnormalised relative abundance and clr-transformed relative abundance.

test_thresholdR
This script was used to test whether there is a performance difference when using different thresholds(0.3,0.5,0.7) for the network construction. SPIEC-EASI and SparCC were performed with all three thresholds for absolute and relative abundance.
