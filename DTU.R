# Differential transcript usage

library(tximport)
library(GenomicFeatures)
library(topGO)
library(DRIMSeq)
library(AnnotationDbi)
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

## setting up the files
gene_GO <- readMappings("TAIR_to_GO.delim")

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

txi = tximport(files, type = "salmon",txOut = T, countsFromAbundance = "scaledTPM")

cts = txi$counts
cts = cts[rowSums(cts) > 0, ] # throw out any transcript with no reads

txdb.filename = "./arabidopsisReference/araport11.201606.annotation.sqlite"
#Make TxDb database (only needed once)
gtf = "./arabidopsisReference/Arabidopsis_thaliana.TAIR10.49.gff3"
txdb = makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)

# Load txDB file
txdb <- loadDb(txdb.filename)

# Make sure dplyr is un
txdf <- AnnotationDbi::select(txdb, keys(txdb, "CDSNAME"), "GENEID", "CDSNAME")
txdf[ is.na(txdf$GENEID), 2] = unlist(strsplit(txdf[ is.na(txdf$GENEID), 1], split = "[.][0-9]*")) # To capture those genes which dont have a geneID/gene name
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

# Range of # of reads across files (in millions)
range(colSums(cts)/1e6)

# There are still 35 genes which appear in one of the data sets but not the other so I am just going to remove them (Maybe ensembl is out of date I have no other idea why this is happening)
cts = cts[rownames(cts) %in% txdf$CDSNAME, ]

txdf <- txdf[match(rownames(cts),txdf$CDSNAME),]
all(rownames(cts) == txdf$TXNAME)

# Set-up for drimSEQ
counts <- data.frame(gene_id=txdf$GENEID,
                     feature_id=txdf$CDSNAME,
                     cts)
# Make sample dataframe
samps = as.data.frame(cbind(rownames(design_full), paste(as.character(design_full[ ,3]), as.character(design_full[, 4]), as.character(design_full[, 5]) ,sep = "_")))
class(samps)

colnames(samps) = c("sample_id", "group")

d <- dmDSdata(counts=counts, samples=samps)
methods(class =class(d))

# set up filtering for DTU, n is the total sample number, n.small is the sample size of the smallest group
n = 36
n.small = 3

d <- dmFilter(d,
min_samps_feature_expr=n.small, min_feature_expr=10, #transcript must have a count of at least 10 in n.small number of samples
min_samps_feature_prop=n.small, min_feature_prop=0.05, #this transcript must have a relative abundance of 5% in at leasr n.small samples
min_samps_gene_expr=n, min_gene_expr=10) # The total count of the corresponding gene must be 10 in all 10 samples

# went from 43,527 transcripts to 16,648 transcripts (with multiple transcripts), genes were only one transcript is observed are removed from this object, which makes sense, since we are testing for differential transcript usage

table(table(counts(d)$gene_id))

factor()
design <- model.matrix(~0+group, data=DRIMSeq::samples(d))
colnames(design)

tmp = d[1:500, ]
16648/500
system.time({
tmp <- dmPrecision(tmp, design=design)
tmp <- dmFit(tmp, design=design)
tmp <- dmTest(tmp, coef="groupm_pst_12")
})
