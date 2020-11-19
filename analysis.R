##Thingy file
library(limma)
library(affy)
##Things to do:
#splitFind needs to be better at determining if it found a real gene (Do ORFs count -> no since ORF's don't appear in desvec) -> the gene could be any of the loci

# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sat Aug 8 18:48:30 EDT 2020



library(Biobase)
library(GEOquery)
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Investigations\\PEN3")
# load series and platform data from GEO

gset <- getGEO("GSE3220", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL198", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
plotMDS(gset@experimentData@name)
# set parameters and draw the plot

dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE3220", '/', annotation(gset), " selected samples", sep ='')
boxplot(exprs(gset), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2) #test if expression is log transformed

#if not log transformed DO THIS!!! If you don't your results will be crazy since limma assumes you are feeding log-transformed data into it so it will output the linear form == insanity
data = exprs(gset) = log2(exprs(gset))

#Function to convert a microarray element to an accession number -> tag the accession number(s) with description(s) and object symbol(s)
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

##Criteria I am aiming for -> looking for genes which are the same between Col-0 Untreated and pen3 untreated qval(0.1)
#pen3_bgh>col_bgh and pen3_ec>col_ec
#qval = 0.05

# #Load the data
# setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Investigations\\PEN3")
# data <- read.delim("pen3 vs wt microarray-modified.txt", comment.char = "!", header = T, as.is = TRUE) 
# colnames(data)[1] = "accession"
# data = data[3:nrow(data), ]
gene_accessions = read.delim("affy_ATH1_array_elements-2010-12-20.txt", comment.char = "!", header = T, as.is = TRUE) 
# ##col-0 Untreated, Col-0 Bgh, Col-0 Ec, Pen3_un, Pen3_Bgh, Pen3_ec
# sub_treat = c(1,9:12,2:4,24, 5:8, 20:23, 13:16, 17:19,25)
# data = data[,sub_treat] #rearrange in my preffered order

#experimental design
genotype = c(rep("col", 11), rep("pen3", 11), "col", "pen3")
length(genotype)
genotype = factor(genotype, levels = c("col", "pen3"))

infection = c(rep("bgh",3), rep("ec", 4), rep("un", 4), rep("bgh", 4), rep("ec", 3), rep("un", 4), "bgh", "ec")
length(infection)
infection = factor(infection, levels = c("un", "bgh", "ec"))

geninf = paste(as.character(genotype), as.character(infection), sep = "_")
length(geninf)
geninf = factor(geninf, levels = c("col_un", "col_bgh", "col_ec", "pen3_un", "pen3_bgh", "pen3_ec"))

model_full <- model.matrix(~0 + geninf) #just use geninf
colnames(model_full) = c("col_un", "col_bgh", "col_ec", "pen3_un", "pen3_bgh", "pen3_ec")


#all the comparisons
fit = lmFit(gset, model_full)

#contrast.matrix <- makeContrasts(pen3_un-col_un, pen3_bgh-col_bgh, pen3_ec-col_ec, levels=model_full)

#contrast.matrix = makeContrasts(col_bgh-col_un, col_ec-col_un, pen3_un-col_un, pen3_bgh-col_un, pen3_ec-col_un, levels = model_full)

contrast.matrix = makeContrasts(pen3_un-col_un, col_ec-col_un, pen3_ec-col_ec, levels = model_full)

contrast.matrix = makeContrasts(col_bgh-col_un, pen3_bgh-pen3_un, pen3_bgh-col_bgh, levels = model_full)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results = decideTests(fit2, adjust.method = "BH", p.value = 0.05)
vennDiagram(results)
vennDiagram(results, include=c("up","down"))

comp = 3
temp = t(as.data.frame(lapply(rownames(results@.Data), splitFind)))
volcanoplot(fit2, coef = comp, highlight= 25, names= temp[,2], main=colnames(fit2$contrasts)[comp])


results@.Data[results@.Data[,1] ==1 & results@.Data[,2]==1 & results@.Data[,3]==1, ]
tab = topTable(fit2, coef=comp, adjust="BH", number = 800)
tab = topTableF(fit2, coef = comp, adjust = "BH", number = 20)

tab = cbind(t(as.data.frame(lapply(rownames(tab), splitFind))), tab)
rownames(tab) = NULL

temp = t(as.data.frame(lapply(rownames(results@.Data[#results@.Data[,1] == 0 & 
          results@.Data[, 2] ==1 & 
          results@.Data[, 3] ==1, ]), splitFind)))
rownames(temp) = NULL

temp = unique(unlist(strsplit(temp, split = ";", fixed = T)))


#Brute forcing to get all the genes in a single excel file with result data
write.fit(fit2, results = results, file = "temp.csv", sep = ",")
data <- read.delim("temp.csv", comment.char = "!", header = T, as.is = TRUE, sep = ",") 
data = data[rownames(data) %in% rownames(results@.Data[ results@.Data[,3] == 1 & results@.Data[,2] ==1 , ]), ]
temp = lapply(rownames(data), splitFind)
temp = t(as.data.frame(temp))
rownames(temp) = NULL
final = cbind(temp[,1:2], data[,2:ncol(data)], temp[, 3])
write.xlsx(final, "comp.xlsx", row.names = F, col.names = T)

sg_y = unique(unlist(strsplit(final[,1], split = ";", fixed = T)))


hello = unlist(lapply(final[,1],function(x){
     x = unlist(strsplit(x, split = ";", fixed = T))
     y = paste0(temp, collapse = "")
     for (i in 1:length(x)){
          x[i] = grepl(pattern = x[i],  y)
     }
     if (sum(as.logical(x))>0){
          return(T)
     } else{
          return(F)
     }
}))



temp = lapply(rownames(fit@.Data[[1]]), splitFind)
temp = t(as.data.frame(temp))

#Make the columns numeric ->Lazy next time use apply()
for( i in 2:ncol(data)){
     data[,i] = as.numeric(data[,i])
}
#calculate average by treatment grou[]
avr = cbind(data[,1], rowMeans(data[, 2:5]), rowMeans(data[, 6:9]), rowMeans(data[, 10:13]), rowMeans(data[, 14:17]), rowMeans(data[, 18:21]), rowMeans(data[, 22:25]))
colnames(avr) = c("acc", "col_un", "col_bgh", "col_ec", "pen3_un", "pen3_bgh", "pen3_ec")

#calculate the foldchange relative to col_un
lfc = apply(avr, MARGIN = 1, function(x){
     tmp = as.numeric(x[3:7])/as.numeric(x[2])
     tmp = c(x[1], tmp)
     return(tmp)
})  
lfc = t(lfc)
colnames(lfc) = colnames(avr)[c(1,3:7)]
lfc = cbind(temp[, 1], temp[, 2], lfc[, 2:6], temp[, 3])

##Actual math -> need to calculate
lfc = lfc[order(as.numeric(lfc[,6]) + as.numeric(lfc[,7]) - as.numeric(lfc[,3]) - as.numeric(lfc[,4]) - as.numeric(lfc[,5]), decreasing = T), ]
rownames(lfc) = NULL
colnames(lfc) = c("accession", "Object symbol", "Col_Bgh", "Col_Ec", "pen3_Un", "pen3_Bgh", "pen3_Ec", "Gene Description")
write.xlsx(lfc, "temp.xlsx", row.names = F, col.names = T)

df = cbind(lfc, desvec[lst]
)



genes.info = function(lst, hpi, linear = T, mock = F){
     hpi = as.character(hpi)
     
     comp = compare.group(hpi = hpi)
     names(comp) = c("ypst_ymg", "mpst_mmg", "mmg_ymg", "mpst_ypst")
     
     res_y= results(allData,contrast = c("group", paste0("ypst", hpi), paste0("ymg", hpi)), alpha = 0.05, pAdjustMethod="BH")
     temp = as.data.frame(res_y@listData)
     rownames(temp) = res_y@rownames
     res_y = temp

     res_m = results(allData,contrast = c("group", paste0("mpst", hpi), paste0("mmg", hpi)),alpha =0.05, pAdjustMethod = "BH")
     temp = as.data.frame(res_m@listData)
     rownames(temp) = res_m@rownames
     res_m = temp
     
     yfc = res_y$log2FoldChange[rownames(res_y) %in% lst]
     mfc = res_m$log2FoldChange[rownames(res_m) %in% lst]
     
     
     if (linear == T){
          yfc = log2linear(yfc)
          mfc = log2linear(mfc)
     }
     if (mock == T){
          res_mock = results(allData,contrast = c("group", paste0("mmg", hpi), paste0("ymg", hpi)),alpha =0.05, pAdjustMethod = "BH")
          temp = as.data.frame(res_mock@listData)
          rownames(temp) = res_mock@rownames
          res_mock = temp
          mockfc = res_mock$log2FoldChange[rownames(res_mock) %in% lst]
          if (linear ==T){
               mockfc = log2linear(mockfc)
          }
          df = cbind(lst, 
               objectSymbol[lst], 
               countMean[lst], 
               yfc, 
               mockfc, 
               mfc, 
               desvec[lst]
          )
          colnames(df) = c("Accession", "Gene name", "Mean count", "Y.Pst-Y.Mock LFC", "M.Mock-Y.Mock LFC", "M.Pst-M.Mock LFC", "Gene description")
          return(df)
     }
     
     df = cbind(lst, 
          objectSymbol[lst], 
          countMean[lst], 
          yfc, 
          res_y$padj[rownames(res_y) %in% lst], 
          mfc, 
          res_m$padj[rownames(res_m) %in% lst], 
          desvec[lst]
     )
     colnames(df) = c("Accession", "Gene name", "Mean count", "Y.Pst-Y.Mock LFC", "q-val", "M.Pst-M.Mock LFC", "q-val", "Gene description")
     return(df)
}