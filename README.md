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

Usage: Usage: CscoreTool-M <windows_Rabl.bed> <input.pair> <OutputPrefix> <N_subcompartments> <N_Rablwindow> <session> [Blacklist.bed] [ExcludedInteractions.txt]
Input parameters

a. windows.bed This file is to specify the genomic windows to analyze. It should be equal-length bed files covering the region of interest, presumably a whole chromosome or whole genome. An example hg19_chr1_10k.bed can be downloaded.These files can also be generatd using the generateEqualLengthBed.cpp program provided here. The chromosome size files used for generating windows.bed are also available for download.

b. input.summary This file is the main input file for Hi-C interactions. We accept the same format as the HiCsummary file format for HOMER runHiCpca.pl. See http://homer.ucsd.edu/homer/interactions/HiCtagDirectory.html An example file test.summary.gz can be downloaded. This is 0.5% randomly selected reads in chr1 from the High-resolution GM12878 cell Hi-C dataset (Rao, 2014).

c. OutputPrefix This is the prefix for output files.

d. session This the number of sessions to use. The number of choice depends on the resource available. 

e. Blacklist.bed This is the encode blacklist file which lists the blacklist regions excluded from ChIP-seq analysis, possibly due to errors or uncertainties in genome assembly or copy number variations between genome. If the program takes 7 arguments, it will take the 7th argument as Blacklisted.bed.

f. ExcludedInteractions.txt This is 

Example run:

CscoreTool1.1 hg19_chr1_10k.bed Test.summary Test_chr1_10k 12 1000000

Output
There are 3 output files.

XXX_cscore.txt This is the Cscore estimated for each genomic window. 
XXX_cscore.bedgraph This is the bedgraph file made for visualization. Ths c-scores are normalized to add up to 1 for each genomic region. Regions with no reads are not shown.
XXX_rabl.txt estimated Rabl effect, shown as a matrix.
