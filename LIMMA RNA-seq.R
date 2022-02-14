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
txdb.filename = "./arabidopsisReference/Arabidopsis_thaliana.TAIR10.52.annotation.sqlite"
#Make TxDb database (only needed once)

# gtf = "./arabidopsisReference/Arabidopsis_thaliana.TAIR10.52.gff3"
# txdb = makeTxDbFromGFF(gtf, format = "gff3", organism = "Arabidopsis thaliana", taxonomyId = 3702)
# saveDb(txdb, txdb.filename)


# Load txDB file
txdb <- loadDb(txdb.filename)

txdf <- AnnotationDbi::select(txdb, keys(txdb, "CDSNAME"), "GENEID", "CDSNAME")
txi = summarizeToGene(txi, tx2gene = txdf, countsFromAbundance = "no") #summarize transcript level data to gene-level
cts = txi$counts
cts = cts[rowSums(cts) > 360, ] # throw out any transcript with no reads
# cts = cts[,c(4:9, 13:18, 22:27, 31:36)]
# design_full = design_full[c(4:9, 13:18, 22:27, 31:36), ]

# Make sample dataframe
samps = as.data.frame(cbind(rownames(design_full), paste(as.character(design_full[ ,3]), as.character(design_full[, 4]), as.character(design_full[, 5]) ,sep = "_")))
class(samps)

colnames(samps) = c("sample_id", "group")




#Transition into LIMMA
dge= DGEList(counts = cts, samples = samps)
design = cbind(factor(design_full[, 3], levels = c("y", "m")), factor(design_full[, 4], levels = c("mg", "pst")), factor(design_full[, 5]))


group = interaction(age,infection,hpi)
mm = model.matrix(~0+group)
colnames(mm) = c("y.mg.0", "m.mg.0", "y.pst.0", "m.pst.0", "y.mg.12", "m.mg.12", "y.pst.12", "m.pst.12", "y.mg.24", "m.mg.24", "y.pst.24", "m.pst.24")
# keep = filterByExpr(dge, design, min.total.count = 360)
# dge = dge[keep,,keep.lib.sizes=F]
# dge = calcNormFactors(dge)

##The samples are relatively close in library size so I'm going to use limma-trend
voomcts = voom(dge, design, plot = T)
voomcts = voom(dge, mm, plot = T)
fit = lmFit(voomcts, mm)
head(coef(fit))

ct.m = makeContrasts(Dif0hr =(m.pst.0-m.mg.0)-(y.pst.0-y.mg.0),
  Dif12hr =(m.pst.12-m.mg.12)-(y.pst.12-y.mg.12), 
  Dif24hr =(m.pst.24-m.mg.24)-(y.pst.24-y.mg.24),
  levels = mm)
fit2 = contrasts.fit(fit, ct.m)
fit2 <- eBayes(fit2)
top.table = topTable(fit2, n = Inf)


objectSymbol[rownames(head(top.table, 80))]

plot