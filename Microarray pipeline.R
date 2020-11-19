##Thingy file
options(java.parameters = "-Xmx2048m")
library(limma)
library(affy)
library(affyPLM)
library(ggplot2)
library(xlsx)
##Things to do:
#splitFind needs to be better at determining if it found a real gene (Do ORFs count -> no since ORF's don't appear in desvec) -> the gene could be any of the loci

# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sat Aug 8 18:48:30 EDT 2020


##Criteria I am aiming for -> looking for genes which are the same between Col-0 Untreated and pen3 untreated qval(0.1)
#pen3_bgh>col_bgh and pen3_ec>col_ec
#qval = 0.05

splitFind = function(x){
        acc = gene_accessions$locus[grep(x, gene_accessions[,1])]
        des = paste(desvec[unlist(strsplit(acc, split = ";"))], collapse = ";")
        sym = objectSymbol[unlist(strsplit(acc, split = ";"))]
        
        sym = paste(sym[is.na(sym)==F], collapse = ";")
        if (sym == ""){
                sym = gene_accessions$locus[grep(x,gene_accessions[,1])]
        }
        if (is.na(des) == T){
                sym = ""
        }
        tmp = c(acc, sym, des)
        return(tmp)
}

#Load the data
setwd("C:\\Users\\Garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop\\2009 microarray")
data = read.xlsx2("ARR-microarray.xlsx",sheetIndex = 1, header = T, startRow = 2) #I believe this is unfiltered (with present call run) normalized data
data = data[ , c(1,2,4,6,8,10,12,3,5,7,9,11,13)]
data[ , 2:7] = as.data.frame(sapply(data[, 2:7], as.numeric))
rownames(data) = data[,1]
data = data[,2:13]
colnames(data) = c("bot0004","bot0005", "bot0006", "bot0007", "bot0008", "bot0009","bot0004_sig","bot0005_sig", "bot0006_sig", "bot0007_sig", "bot0008_sig", "bot0009_sig")

#Remove genes which are not present in either all Mock samples or all Pst-treated samples

for (i in 7:ncol(data)){
        data[, i] = data[, i] == "P"
} #Any samples labelled P are labelled true


oldGOI = rowSums(data[,7:12]) >=3  #cameron 2009 paper-> all 'good' probes not A in any sample
genesOI = rowSums(data[,c(7,9,10)]) > 1 & rowSums(data[,c(8,11,12)]) > 1 #based on what I think
# genesOI = rowSums(data[,c(7,9,10)]) == 3 | rowSums(data[,c(8,11,12)]) ==3
data = data[genesOI, 1:6]


# #remove affy data
# data = data[65:nrow(data),]

#Phenotype data (rows are samples, columns are covariates)
pData = as.data.frame(cbind(c("mock", "pst", "mock", "mock", "pst", "pst"),c(T,T,T,F,F,T)))
# pData = as.data.frame(as.factor(c("mock", "pst", "mock", "pst", "mock", "pst")))
rownames(pData) = colnames(data)
colnames(pData) = c("treatment", "known")
metadata = data.frame(labelDescription = c("Treatment", "Do we know the treatment type"), row.names = c("treatmennt","known") )
pData = new("AnnotatedDataFrame", data = pData, varMetadata = metadata)

#Experiment data
eData = new("MIAME", name = "Garrett Nunn", lab = "Robin Cameron Lab", contact = "rcamero@mcmaster.ca", title = "ARR Microarray Experiment")

#import all data into an expressionSet object
data = as.matrix(data)
minSet = ExpressionSet(assayData = data, phenoData = pData, experimentData = eData)

sampleNames(minSet)

#Load gene_accessions
gene_accessions = read.delim("affy_ATH1_array_elements-2010-12-20.txt", comment.char = "!", header = T, as.is = TRUE) 

#normalization
# eset = normalize.ExpressionSet.quantiles(minSet, transfn = "none")
eset = normalize.ExpressionSet.loess(minSet, transfn = "none")# this one


#Pre-processing

tiff(filename = "microarraynormalisedExpression.tiff", height = 1000, width = 1000)
par(mar=c(2+round(max(nchar(colnames(data)))/2),4,2,1))
boxplot(exprs(eset), boxwex=0.7, notch=T, outline=FALSE, las=2)
boxplot(exprs(minSet), boxwex=0.7, notch=T, outline=FALSE, las=2)
dev.off()

pca <- prcomp(t(exprs(eset)))
percentVar <- pca$sdev^2/sum(pca$sdev^2)


#Selecting the principle components
pc_x = 1
pc_y = 2
d <- data.frame(PC1 = pca$x[, pc_x], PC2 = pca$x[, pc_y], 
                group = colnames(data),
                test = c("mock", "pst", "mock", "mock", "pst", "pst")
)
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + geom_point(size = 4, stroke = 1.5) +  xlab(paste0("PC ",pc_x," (", round(percentVar[pc_x] * 100), "% variance)")) + ylab(paste0("PC ",pc_y," (", round(percentVar[pc_y] * 100), "% variance)")) + 
        coord_fixed() +theme(panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank(), 
                             axis.line = element_line(colour = "black", size=1),
                             axis.title.x=element_text(size=15),
                             #axis.text.x=element_blank()),
                             axis.ticks=element_line(colour = "black", size =1),
                             axis.ticks.length = unit(5,"points") ,
                             axis.title.y = element_text(size=15),
                             legend.position = "right",
                             axis.text = element_text(size=15),
        )

tiff(filename = "PCA_microarray.tiff", height = 2000, width = 2000) #opens a tiff device
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + geom_point(size = 20, stroke = 6) +  xlab("") + ylab("") + coord_fixed() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size=4), axis.title.x=element_text(size=15), #axis.text.x=element_blank()),
                                                                                                                                                      axis.ticks=element_line(colour = "black", size =4),
                                                                                                                                                      axis.ticks.length = unit(20,"points") ,
                                                                                                                                                      axis.title.y = element_text(size=15),
                                                                                                                                                      legend.position = "right", axis.text = element_text(size= 75)
)
dev.off()



treatment = eset@phenoData@data$treatment

model_full <- model.matrix(~0 + treatment)
colnames(model_full) = c("mock", "pst")

fit = lmFit(log2(exprs(eset)), model_full) #requires log expression values

contrast.matrix <- makeContrasts(pst-mock, levels=model_full)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
temp = treat(fit, lfc = 1)
        
results = decideTests(fit2, adjust.method = "none",lfc = log2(1.5), p.value = 0.15)
vennDiagram(results)

genes = lapply( rownames(fit$coefficients),splitFind)
genes = t(as.data.frame(genes))

volcanoplot(fit2, highlight = 10, names= genes[,2], main="pst-mock")

tab = topTable(fit2, adjust.method = "none",lfc = log2(1.5), p.value = 0.05, resort.by = "logFC", number = 1000)
tab = cbind(tab, t(as.data.frame(lapply(rownames(tab), splitFind)))[,2:1])




#Compare to RNA-seq data
microGenes = tab[as.numeric(tab$logFC) > 0,8] #Microarray mature Pst induced genes
#Split up accession numbers
microGenes = unique(unlist(strsplit(microGenes, split = ";")))


tmp = compare.group(hpi = "12") #collect RNA-seq data

#Young Pst-induced genes
young = tmp$ypst12ymg12
young = young[young$log2FoldChange > log2(1.5) & young$padj < 0.05, ]
young = rownames(young)

#Mature Pst-induced genes
mature = tmp$mpst12mmg12
mature = mature[mature$log2FoldChange > log2(1.5) & mature$pvalue < 0.05, ]
mature = rownames(mature)

mature = find.unique(hpi = "12")
