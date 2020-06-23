##Original gene_association, Gene-to-GO, and ObjectSymbol set-up
library(DESeq2)
library(RCurl)
library(tximport)
library(readr)
library("RColorBrewer")
library(gplots)
library(topGO) 
library(Rgraphviz) 
library(dplyr)
library(ggplot2)
library(colorspace)
library(ggrepel)
library(VennDiagram)

gene_associations <- read.delim("gene_association_final.txt", comment.char = "!", header = FALSE, as.is = TRUE) 
colnames(gene_associations) <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
                                 "DB:Reference", "Evidence", "With_From", "Aspect", "DB_Object_Name",
                                 "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_by") 
#I didn't want to look at locus ID's from TAIR so I took the first item of DB_Object_Synonym which is the gene ID and put it in the DB_Object_ID column as that is more useful
gene_associations$DB_Object_ID = res
gen_association_save <- sapply(gene_associations, function(x){paste(x, collapse = ", ")}) 
gen_association_save <- data.frame("Gene Association" = gene_associations)
write.table(gen_association_save, file = "gene_association_final.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = F)

# Trimming the dataframe so that it's only what we're interested in
gene_associations <- gene_associations[,c(2,3,5,7,9,14)]



# Go through every unique gene and pull out any GO_ID associated with the given gene then give
gene_GO <- lapply(unique(gene_associations$DB_Object_ID), function(x){tmp <- gene_associations %>% filter(DB_Object_ID == x)
return(tmp$GO_ID)})
names(gene_GO) <- unique(gene_associations$DB_Object_ID) 
#Save it
gene_GO_save <- sapply(gene_GO, function(x){paste(x, collapse = ", ")}) 
gene_GO_save <- as.data.frame( gene_GO_save) # Making a dataframe
write.table(gene_GO_save, file = "TAIR_to_GO.delim", sep = "\t", quote = FALSE, col.names = FALSE) 

objectSymbol = lapply(unique(gene_associations$DB_Object_ID), function(x){tmp <- gene_associations %>% filter(DB_Object_ID == x)
return(tmp$DB_Object_Symbol)}) ##If gene assoc
names(objectSymbol) = unique(gene_associations$DB_Object_ID)
objectSymbol = unlist2(objectSymbol)
objectSymbol = objectSymbol[unique(names(objectSymbol))]

#load Pre-meade Gene-to-go from file
gene_GO <- readMappings("TAIR_to_GO.delim")

setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")
files <- file.path("counts", list.files("counts"))


samples <- c("m_mg_0h_s1", "m_mg_0h_s2","m_mg_0h_s3","m_mg_12h_s1","m_mg_12h_s2","m_mg_12h_s3","m_mg_24h_s1","m_mg_24h_s2","m_mg_24h_s3","m_pst_0h_s1","m_pst_0h_s2","m_pst_0h_s3","m_pst_12h_s1","m_pst_12h_s2","m_pst_12h_s3","m_pst_24h_s1","m_pst_24h_s2","m_pst_24h_s3","y_mg_0h_s1","y_mg_0h_s2","y_mg_0h_s3","y_mg_12h_s1","y_mg_12h_s2","y_mg_12h_s3","y_mg_24h_s1","y_mg_24h_s2","y_mg_24h_s3","y_pst_0h_s1","y_pst_0h_s2","y_pst_0h_s3","y_pst_12h_s1","y_pst_12h_s2","y_pst_12h_s3","y_pst_24h_s1","y_pst_24h_s2","y_pst_24h_s3")
names(files) <- samples

vg = function(accession){
  plotCounts(all.data2, 
             gene = toupper(accession),
             intgroup="group",
             pch = 20, col = "red")
}