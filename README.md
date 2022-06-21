# CscoreTool-M
Multiple sub-compartment analysis from Hi-C data

Purpose

This program is to do multiple sub-compartment analysis for Hi-C data. Test data are also provided. Several other files are also available for download, including the bedgraph files involved in the CscoreTool-M paper.

Installation

The program is simple one-file executable for linux system. Just download CscoreTool-M Then run the following command:

chmod +x CscoreTool-M

Then you can type ./CscoreTool-M to run it.

The executable is compiled on a CentOS linux machine. If it doesn't work, download the CscoreTool-M.cpp, twister.h and twister.cpp files and put them in the same folder. Then run

g++ CscoreTool-M.cpp twister.cpp -fopenmp -O3 -o CscoreTool-M

You'll get an execuable file CscoreTool-M.

Usage: Usage: CscoreTool-M <windows_Rabl.bed> <input.summary> <OutputPrefix> <N_subcompartments> <N_Rablwindow> <session> [Blacklist.bed] [ExcludedInteractions.txt]

Input parameters

a. windows.bed This file is to specify the genomic windows to analyze. It should be an equal-length bed file with the fourth column being the Rabl position, i.e. the relative position between centromere (-0.5) and telomere (0.5). 

b. input.summary This file is the main input file for Hi-C interactions. We accept the same format as the HiCsummary file format for HOMER runHiCpca.pl. See http://homer.ucsd.edu/homer/interactions/HiCtagDirectory.html An example file test.summary.gz can be downloaded. This is 0.5% randomly selected reads in chr1 from the High-resolution GM12878 cell Hi-C dataset (Rao, 2014).

c. OutputPrefix This is the prefix for output files.

d. session This the number of sessions to use. The number of choice depends on the resource available. 

e. Blacklist.bed This is the encode blacklist file which lists the blacklist regions excluded from ChIP-seq analysis, possibly due to copy number variations. If the program takes 7 arguments, it will take the 7th argument as Blacklisted.bed.

f. ExcludedInteractions.txt This is the file for excluding regions having translocations. The format is like
  
chrName1a_chrStart1a_chrEnd1a chrName1b_chrStart1b_chrEnd1b
  
chrName2a_chrStart2a_chrEnd2a chrName2b_chrStart2b_chrEnd2b chrName2c_chrStart2c_chrEnd2c
  
  ......
  
For each line, several genomic regions are listed in the chrName_chrStart_chrEnd format and separated by tab. Then any interactions between genomic regions in the same line will be excluded from the analysis. Note that the expectation are also excluded in the calculation, so it's different from simply delete the interactions from the input file. If the program takes 8 arguments, it will take the 8th argument as ExcludedInteractions.txt.

Example run:

./CscoreTool-M hg38_100k_Rabl.bed Test.summary Test_100k_1 5 20 12 

./CscoreTool-M hg38_100k_Rabl.bed Test.summary Test_100k_2 5 20 12 hg38_blacklisted.bed

./CscoreTool-M hg38_100k_Rabl.bed Test.summary Test_100k_3 5 20 12 hg38_blacklisted.bed Test_excluded.txt

Output
There are 3 output files.

XXX_cscore.txt This is the Cscore estimated for each genomic window. 
XXX_cscore.bedgraph This is the bedgraph file made for visualization. Ths c-scores are normalized to add up to 1 for each genomic region. Regions with no reads are not shown.
XXX_rabl.txt estimated Rabl effect, shown as a matrix.
