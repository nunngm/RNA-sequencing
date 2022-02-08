#LIMMA for ARR data
library(edgeR)
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(dplyr)
set.seed(31138)

#laptop directory
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")

#Desktop directory
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")


gene_associations <- read.delim("gene_association_final.txt", comment.char = "!", header = FALSE, as.is = TRUE) 
colnames(gene_associations) <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
                                 "DB:Reference", "Evidence", "With_From", "Aspect", "DB_Object_Name",
                                 "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_by")

gene_associations <- gene_associations[,c(2,3,5,7,9,14)]


#Collect filenames and label samples
files = file.path(list.dirs("salmon_quant", recursive=F)[], "quant.sf") #[-1] is to ignore the salmon_quant directory

# Samples in alphabetical order
samples <- c("m_mg_0h_s1", "m_mg_0h_s2","m_mg_0h_s3","m_mg_12h_s1","m_mg_12h_s2","m_mg_12h_s3","m_mg_24h_s1","m_mg_24h_s2","m_mg_24h_s3","m_pst_0h_s1","m_pst_0h_s2","m_pst_0h_s3","m_pst_12h_s1","m_pst_12h_s2","m_pst_12h_s3","m_pst_24h_s1","m_pst_24h_s2","m_pst_24h_s3","y_mg_0h_s1","y_mg_0h_s2","y_mg_0h_s3","y_mg_12h_s1","y_mg_12h_s2","y_mg_12h_s3","y_mg_24h_s1","y_mg_24h_s2","y_mg_24h_s3","y_pst_0h_s1","y_pst_0h_s2","y_pst_0h_s3","y_pst_12h_s1","y_pst_12h_s2","y_pst_12h_s3","y_pst_24h_s1","y_pst_24h_s2","y_pst_24h_s3")
names(files) <- samples

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


#tximport
txi = tximport(files, type = "salmon",txOut = T, countsFromAbundance = "no")

## Generate tx2gene file
txdb.filename = "./arabidopsisReference/araport11.201606.Jan2022.annotation.sqlite"
#Make TxDb database (only needed once)
gtf = "./arabidopsisReference/Araport11_GFF3_genes_transposons.Jan032022.gff"
txdb = makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)

# Load txDB file
txdb <- loadDb(txdb.filename)

txi = summarizeToGene(txi, tx2gene = txdf, countsFromAbundance = "no") #summarize transcript level data to gene-level
cts = txi$counts
cts = cts[rowSums(cts) > 0, ] # throw out any transcript with no reads

# Make sample dataframe
samps = as.data.frame(cbind(rownames(design_full), paste(as.character(design_full[ ,3]), as.character(design_full[, 4]), as.character(design_full[, 5]) ,sep = "_")))
class(samps)

colnames(samps) = c("sample_id", "group")


##Very general model that just has all the sample information
model_full <- formula(~age+infection+hpi)
rawData <- DESeqDataSetFromHTSeqCount(design_full,design=model_full)
counts = counts(rawData)

#Transition into LIMMA
dge= DGEList(counts = cts, samples = samps)
design = cbind(factor(design_full[, 3], levels = c("y", "m")), factor(design_full[, 4], levels = c("mg", "pst")), factor(design_full[, 5]))
keep = filterByExpr(dge, design, min.total.count = 360)
dge = dge[keep,,keep.lib.sizes=F]
dge = calcNormFactors(dge)

##The samples are relatively close in library size so I'm going to use limma-trend
logCPM = cpm(dge, log = T, prior.count = 3)
fit = lmfit(logCPM, design)