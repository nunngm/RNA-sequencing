#LIMMA for ARR data
library(edgeR)
library(DESeq2)
set.seed(31138)

#laptop directory
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")

#Desktop directory
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")

files <- file.path("counts", list.files("counts"))
samples <- c("m_mg_0h_s1", "m_mg_0h_s2","m_mg_0h_s3","m_mg_12h_s1","m_mg_12h_s2","m_mg_12h_s3","m_mg_24h_s1","m_mg_24h_s2","m_mg_24h_s3","m_pst_0h_s1","m_pst_0h_s2","m_pst_0h_s3","m_pst_12h_s1","m_pst_12h_s2","m_pst_12h_s3","m_pst_24h_s1","m_pst_24h_s2","m_pst_24h_s3","y_mg_0h_s1","y_mg_0h_s2","y_mg_0h_s3","y_mg_12h_s1","y_mg_12h_s2","y_mg_12h_s3","y_mg_24h_s1","y_mg_24h_s2","y_mg_24h_s3","y_pst_0h_s1","y_pst_0h_s2","y_pst_0h_s3","y_pst_12h_s1","y_pst_12h_s2","y_pst_12h_s3","y_pst_24h_s1","y_pst_24h_s2","y_pst_24h_s3")
names(files) <- samples

setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")

#Desktop directory
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")


gene_associations <- read.delim("gene_association_final.txt", comment.char = "!", header = FALSE, as.is = TRUE) 
colnames(gene_associations) <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
                                 "DB:Reference", "Evidence", "With_From", "Aspect", "DB_Object_Name",
                                 "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_by")

# # Uncomment if this is the first time running this code this helps produce a gene_association file to work with in the future
# #I didn't want to look at locus ID's from TAIR so I took the first item of DB_Object_Synonym which is the gene ID and put it in the DB_Object_ID column as that is more useful
# gene_associations$DB_Object_ID = res
# gen_association_save <- sapply(gene_associations, function(x){paste(x, collapse = ", ")}) 
# gen_association_save <- data.frame("Gene Association" = gene_associations)
# write.table(gen_association_save, file = "gene_association_final.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = F)

# Trimming the dataframe so that it's only what we're interested in
gene_associations <- gene_associations[,c(2,3,5,7,9,14)]


##Make Gene-Go object from scratch
# # Go through every unique gene and pull out any GO_ID associated with the given gene then give
# gene_GO <- lapply(unique(gene_associations$DB_Object_ID), function(x){tmp <- gene_associations %>% filter(DB_Object_ID == x)
# return(tmp$GO_ID)})
# names(gene_GO) <- unique(gene_associations$DB_Object_ID) 
# #Save it
# gene_GO_save <- sapply(gene_GO, function(x){paste(x, collapse = ", ")}) 
# gene_GO_save <- as.data.frame( gene_GO_save) # Making a dataframe
# write.table(gene_GO_save, file = "TAIR_to_GO.delim", sep = "\t", quote = FALSE, col.names = FALSE) 


#Load a pre-made Gene-to-go from file
gene_GO <- readMappings("TAIR_to_GO.delim")

#Collect filenames and label samples
files <- file.path("counts", list.files("counts"))
samples <- c("m_mg_0h_s1", "m_mg_0h_s2","m_mg_0h_s3","m_mg_12h_s1","m_mg_12h_s2","m_mg_12h_s3","m_mg_24h_s1","m_mg_24h_s2","m_mg_24h_s3","m_pst_0h_s1","m_pst_0h_s2","m_pst_0h_s3","m_pst_12h_s1","m_pst_12h_s2","m_pst_12h_s3","m_pst_24h_s1","m_pst_24h_s2","m_pst_24h_s3","y_mg_0h_s1","y_mg_0h_s2","y_mg_0h_s3","y_mg_12h_s1","y_mg_12h_s2","y_mg_12h_s3","y_mg_24h_s1","y_mg_24h_s2","y_mg_24h_s3","y_pst_0h_s1","y_pst_0h_s2","y_pst_0h_s3","y_pst_12h_s1","y_pst_12h_s2","y_pst_12h_s3","y_pst_24h_s1","y_pst_24h_s2","y_pst_24h_s3")
names(files) <- samples

#Subset dataset by time
# sub_0 = c(1:3,10:12,19:21,28:30)
# sub_12 = c(4:6,13:15,22:24,31:33)
# sub_24 = c(7:9,16:18,25:27,34:36)

#ymg ypst mmg mpst
sub_0 = c(19:21,28:30, 1:3,10:12) 
sub_12 = c(22:24,31:33,4:6,13:15)
sub_24 = c(25:27,34:36, 7:9,16:18)
##Sample information to build a model
infection = factor(c(rep("mg", 9), rep("pst", 9), rep("mg", 9), rep("pst", 9)), levels = c("mg", "pst")) #Treatment type, mg = mock solution, pst = P syringae pv tomato
hpi = factor(rep(c(rep(0, 3), rep(12, 3), rep(24, 3)), 4), levels = c(0, 12, 24)) #Hours post infiltration
age = factor(c(rep("m",18),rep("y",18)), levels = c("y", "m")) #Plant age y = 3.5 wpg, m = 6.5 wpg

design_full <- data.frame(sample=names(files),
                          file=files,
                          age=age,
                          infection =infection,
                          hpi=hpi
                )
design_full

##Very general model that just has all the sample information
model_full <- formula(~age+infection+hpi)
rawData <- DESeqDataSetFromHTSeqCount(design_full,design=model_full)
counts = counts(rawData)

#Transition into LIMMA
dge= DGEList(counts = counts)
design = cbind(factor(design_full[, 3], levels = c("y", "m")), factor(design_full[, 4], levels = c("mg", "pst")), factor(design_full[, 5]))
keep = filterByExpr(dge, design, min.total.count = 360)
dge = dge[keep,,keep.lib.sizes=F]
dge = calcNormFactors(dge)

##The samples are relatively close in library size so I'm going to use limma-trend
logCPM = cpm(dge, log = T, prior.count = 3)
fit = lmfit(logCPM, design)