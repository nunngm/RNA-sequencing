# cisElementEnrichmemt
### The goal of this is to take the output of "Correlation matrix.R" and use a specific co-expression network and then take the upstream gene FASTA file and output two files. 
# File 1 upstream region of genes of interest
# File 2 upstream region of all other genes
# Gong to start with a gene network which is uniquely up-regulated in mature plants in response to infection
library(Biostrings)
library(dplyr)

setwd("C:\\Users\\Garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-Cis-element enrichment")

data = as.character(readDNAStringSet("TAIR10_upstream_500_20101028.txt"))
names(data) = t(as.data.frame(strsplit(names(data), split = " |", fixed = T)))[, 1]


genesOI = clustcolours[clustcolours %in% "purple"]
controlGenes = clustcolours[! clustcolours %in% "purple"]

genesOI = data[names(data) %in% names(genesOI)]
controlGenes = data[names(data) %in% names(controlGenes)]


writeFasta = function(sequences, filename) {
  lines2write = c()
  for (i in 1:length(sequences)){
    lines2write = c(lines2write, paste(">", names(sequences)[i], sep = " "))
    lines2write = c(lines2write, sequences[ i])
  }
  writeLines(lines2write, filename)
}

setwd("C:\\Users\\Garrett\\Documents\\MobaXterm\\home")

writeFasta(genesOI, "genesOI.txt")
writeFasta(controlGenes, "controlGenes.txt")
