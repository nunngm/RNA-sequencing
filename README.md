# R scripts for the analysis of our RNAseq data
Code to accompany ARR transcriptomics paper

## Information on the files
- `ARR-RNA-seqsetup.R`: Sets up the basic experimental variables and functions for DESeq2 analysis
- `Gene_association Gene2Go and objectSymbol setup.R`: Antiquated. Important parts have been included in other files
- `Heatmaps.R`: Contains code for all heatmaps and downstream analysis
- `Isolating Genes.R`: Most of this code has be incorporated into the find.unique() function in the setup file
- `Microarray pipeline.R`: Pipeline for analyzing the microarray data from the 2009 paper "Forward and reverse genetics to identify genes involved in the age‐related resistance response in Arabidopsis thaliana"
- `PCA.R`: Contains code for all PCAs and accompanying plots
- `The One.Rmd`: Antiquated (was split into the array of files we have now)
- `TimeCourse_analysis.R`: Contains code to identify genes which had differential expression patterns over time in mature plants
- `Venn diagrams and TopGO.R`: Contains code for Venn diagrams and TopGo analysis
- `Volcano plot.R`: Contains code for producing volcano plots and subsequent analysis
- `analysis.R`: An antiquated pipeline for analyzing microarray data
