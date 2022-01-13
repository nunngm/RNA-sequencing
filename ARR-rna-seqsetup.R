##RNA-sequencing analysis for ARR
#As a note before the coding begins, throughout the file I refer to samples which were collected at 0.25 hours-post infiltration (hpi) as 0. This was just to simplify the coding process
options(java.parameters = "-Xmx2048m")
options(java.parameters = "-Xmx1024m")
library(DESeq2)
library(RCurl)
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
library(scales)
library(xlsx)
library(stringi)
set.seed(31138)

#laptop directory
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


## Grouping all the experimental variables into distinct treatment groups
rawData$group <- factor(paste0(rawData$age,rawData$infection,rawData$hpi))
rawData$group = factor(rawData$group,levels=c("ymg0","ymg12","ymg24","ypst0","ypst12","ypst24","mmg0","mmg12","mmg24","mpst0","mpst12","mpst24"))

# Filter lowly expressed genes
keep <- rowMeans(counts(rawData)) >= 10 #Genes which on average have less than 10 reads
#keep = rowMeans(counts(rawData)[, 1:18]) >=10 &rowMeans(counts(rawData)[, 19:36]) >=10
rawData <- rawData[keep,]

# Build DESeq object based on distinct treatment groups
rawData@design = ~group
allData <- DESeq(rawData)

##A really nice function (if a little messy) that completes only the comparisons you deem relevent but does not compare between timepoints
        #ie if you set the age to "y" then it will complete ypst - ymg at every time point
compare.group = function(age = c("y", "m"), infection = c("mg", "pst"), hpi = c("0", "12", "24"), altH = c("greaterAbs", "lessAbs", "greater", "less")){
     compGroups = lapply(age, function(x){tmp = paste0(x,infection)
                    return(tmp)
           })
     for (i in 1:length(compGroups)){
          compGroups[[i]] = lapply(compGroups[[i]], function(x){tmp = paste0(x, hpi)
                              return(tmp)
          })
     }
     out = character()
     for (i in 1:length(compGroups[[1]][[1]])){
          if(length(compGroups[[1]]) > 1){
               for (j in 1:length(compGroups)){
                    out = rbind(out, c(compGroups[[j]][[2]][i], compGroups[[j]][[1]][i]))
               }
          }
          if (length(compGroups) > 1){
               for (j in 1:length(compGroups[[1]])){
                    out = rbind(out, c(compGroups[[2]][[j]][i], compGroups[[1]][[j]][i]))
               }
          }
     }
     ##need to get rownames efficiently
     final = list()
     for (i in 1:nrow(out)){
          tmp = results(allData,contrast = c("group", out[i, 1], out[i, 2]), alpha = 0.05, pAdjustMethod="BH", altHypothesis = altH)
          tmp = as.data.frame(tmp@listData, row.names = tmp@rownames)
          final[[i]] = tmp
     }
     names(final) = paste0(out[,1], out[,2])
     return(final)
}

# Set-up the main biologically relevent comparisons
# contrast = c(intgroup,Var1, Var2) -> contrasts Var 1 to reference (Var2)
# pAdjustMethod ="BH" -> adjust the p-value based on the false discovery rate
##The compare.group() function does a much better job of this -> use that
res_y= results(allData,contrast = c("group", "ypst24", "ymg24"), alpha = 0.05, pAdjustMethod="BH")
res_y = as.data.frame(res_y@listData, row.names = res_y@rownames)

res_m = results(allData,contrast = c("group","mpst24","mmg24"),alpha =0.05, pAdjustMethod = "BH")
res_m = as.data.frame(res_m@listData, row.names = res_m@rownames)

res_mock = results(allData,contrast = c("group","mmg24","ymg24"),alpha =0.05, pAdjustMethod = "BH")
res_mock = as.data.frame(res_mock@listData, row.names = res_mock@rownames)

res_pst = results(allData, contrast = c("group", "mpst24", "ypst24"), alpha = 0.05, pAdjustMethod = "BH")
res_pst = as.data.frame(res_pst@listData, row.names = res_pst@rownames)

## Gather gene descriptions from gene_description file on TAIR
gene_descriptions = read.delim("gene_description.delim",sep = "\t",header = FALSE,stringsAsFactors = F,quote = "") 
gene_descriptions[,1]=tools::file_path_sans_ext(gene_descriptions[,1]) #removes ".1 or.2 ...etc" after each gene (these indicate the splice variant) 
gene_descriptions = gene_descriptions[,c(1,4)] #remove unnecessary columns
colnames(gene_descriptions) = c("accession", "description")
gene_descriptions = gene_descriptions %>% filter(!description == "") #remove splice variants with empty description
gene_descriptions = gene_descriptions %>% distinct(accession, .keep_all = T) #keep only the first splice variant with a description

#Makes a vector similar to objectSymbol to go from accession number to description
desvec = gene_descriptions[, 2]
names(desvec) = gene_descriptions[, 1] 

## objectSymbol is a great object that takes you from accession number to the first object symbol for example EDS16 (aka SID2/ICS1)
objectSymbol = lapply(unique(gene_associations$DB_Object_ID), function(x){tmp <- gene_associations %>% filter(DB_Object_ID == x)
return(tmp$DB_Object_Symbol)}) ##If gene assoc
names(objectSymbol) = unique(gene_associations$DB_Object_ID)
objectSymbol = unlist2(objectSymbol)
objectSymbol = objectSymbol[unique(names(objectSymbol))]

##Wrapper for objectSymbol to use on a list of genes
getGeneName = function(x){
     temp = objectSymbol[x]
     if (is.na(temp)==T){
          return(x)
     }
     return(temp)
}

countMean = res_y$baseMean
names(countMean) = rownames(res_y)

## For a given list of genes collect all the relevent information at a given time point
genes.info = function(lst, hpi, linear = T, mock = F){
     hpi = as.character(hpi)
     
     comp = compare.group(hpi = hpi)
     names(comp) = c("ypst_ymg", "mpst_mmg", "mmg_ymg", "mpst_ypst")
     
     res_y= results(allData,contrast = c("group", paste0("ypst", hpi), paste0("ymg", hpi)), alpha = 0.05, pAdjustMethod="BH")
     res_y = as.data.frame(res_y@listData, row.names = res_y@rownames)

     res_m = results(allData,contrast = c("group", paste0("mpst", hpi), paste0("mmg", hpi)),alpha =0.05, pAdjustMethod = "BH")
     res_m = as.data.frame(res_m@listData, row.names = res_m@rownames)
     
     yfc = res_y$log2FoldChange[rownames(res_y) %in% lst]
     mfc = res_m$log2FoldChange[rownames(res_m) %in% lst]
     
     
     if (linear == T){
          yfc = log2linear(yfc)
          mfc = log2linear(mfc)
     }
     if (mock == T){
          res_mock = results(allData,contrast = c("group", paste0("mmg", hpi), paste0("ymg", hpi)),alpha =0.05, pAdjustMethod = "BH")
          res_mock = as.data.frame(res_mock@listData, row.names = res_mock@rownames)

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
# For a given list of genes collect all the relevent information at a given time point
genes.info = function(lst, hpi, linear = T, mock = F)

## Takes an gene accession number and displays/saves the plot of gene expression over time
view.gene = function(accession, fileName = objectSymbol[toupper(accession)], graph = F){
     geneCount <- plotCounts(allData, gene = toupper(accession), 
                        intgroup = c("age", "infection","hpi"), returnData = TRUE)
     geneCount$hpi <- as.numeric(as.character(geneCount$hpi))
     geneCount$hpi[geneCount$hpi == 0] = 0.25
     geneCount = cbind(geneCount, paste0(as.character(geneCount$age), as.character(geneCount$infection)))
     colnames(geneCount)[ncol(geneCount)] = "ageinf"
     
     p = ggplot(geneCount,
       aes(x = hpi, y = count, color = factor(age, levels = c("y", "m")), lty = factor(infection, levels = c("mg", "pst")), group = ageinf )) +
          #geom_point() + #displays individual counts
          stat_summary(fun=mean, geom="line", size = 1) + 
          stat_summary(fun = mean,
                    fun.min = function(x) {mean(x) - sd(x)}, 
                    fun.max = function(x) {mean(x) + sd(x)}, 
                    geom = "pointrange", lty =1 , size =1)+theme(panel.grid.major = element_blank(), 
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
         legend.text = element_text(size=15)
          ) + ggtitle(fileName)
     if (graph == T){
          ggsave(paste0(fileName, ".png"), plot = p)
          
          p

     } else {
         p  
     }
}

# Simply takes a bunch of log2FoldChanges and returns the linear fold change
log2linear = function(lfc){
     lfc = lapply(lfc, function(x){
          temp = 2^(as.double( x))
          if (temp < 1){
               temp = -1/temp
          }
          return(temp)
     })
     return(unlist(lfc))
}

# Takes a path from the clipboard with '\' characters (the way you copy from windows) and returns R-usable file path
get.path = function(){
     return(paste(scan("clipboard", what = "string", sep = "\\"), collapse = "/"))
}
#s
## Finds uniquely up-(or down) regulated genes
find.unique = function(hpi, reverse = F, down = F ){
     comp = compare.group(hpi = hpi)
     names(comp) = c("ypst_ymg", "mpst_mmg", "mmg_ymg", "mpst_ypst")
     genesOI = cbind( comp$ypst_ymg$log2FoldChange, comp$mmg_ymg$log2FoldChange, comp$mpst_mmg$log2FoldChange)
     pvals = cbind(comp$mmg_ymg$padj, comp$mpst_mmg$padj, comp$ypst_ymg$padj)
     
     # for (i in 1:nrow(pvals)){
     #      pvals[i, 1] = -2*sum(log(pvals[i, ]))
     #      pvals[i,1] = 1-pchisq(pvals[i,1], 2*length(pvals[i, ]))
     # }
     genesOI = cbind(rownames(comp$ypst_ymg), genesOI, pvals[ ,1], pvals[, 2], pvals[, 3])
     colnames(genesOI) = c("accession", "ypst_ymg", "mmg_ymg", "mpst_mmg", "qval", "qmpst", "qypst")
     genesOI = as.data.frame(genesOI)
     for (i in 2:ncol(genesOI)){
          genesOI[, i] = as.numeric(as.character(genesOI[, i]))
     }
     
     if (reverse == F & down == F){
          # genesOI = subset(genesOI, qval < 0.01 | qmpst < 0.01)
          genesOI = subset(genesOI, qmpst <0.1)
          genesOI = subset(genesOI, mpst_mmg > 0)
          genesOI = subset(genesOI, ypst_ymg < 0 | qypst > 0.1)
          genesOI = subset(genesOI, mmg_ymg + mpst_mmg - ypst_ymg > log2(1.5))
          genesOI = genesOI[order(genesOI$mmg_ymg + genesOI$mpst_mmg - genesOI$ypst_ymg, decreasing = T), ]
          genesOI = cbind(genesOI, objectSymbol[genesOI$accession], desvec[genesOI$accession])
     } else if (reverse == T & down == F){
          genesOI = subset(genesOI, qypst < 0.1)
          genesOI = subset(genesOI, ypst_ymg > 0)
          genesOI = subset(genesOI, mpst_mmg < 0 | qmpst >0.1)
          genesOI = subset(genesOI, ypst_ymg - mmg_ymg - mpst_mmg > log2(1.5))
          genesOI = genesOI[order(genesOI$ypst_ymg - genesOI$mmg_ymg - genesOI$mpst_mmg, decreasing = T), ]
          genesOI = cbind(genesOI, objectSymbol[genesOI$accession], desvec[genesOI$accession])
     } else {
             genesOI = subset(genesOI, qmpst < 0.05)
             genesOI = subset(genesOI, mpst_mmg < 0)
             genesOI = subset(genesOI, ypst_ymg > 0 | qypst >0.05)
             genesOI = genesOI[order(genesOI$mmg_ymg + genesOI$mpst_mmg - genesOI$ypst_ymg, decreasing = F), ]
             genesOI = cbind(genesOI, objectSymbol[genesOI$accession], desvec[genesOI$accession])
     }

     return(genesOI)
}

## Collect data for volcano plot 
find.volcano = function(hour, n = 15, de = 0, lfc.cut = log2(1.5), adjp.cut = 0.05, SUGs = F) {
        temp = compare.group(hpi= hour)
        data = temp[[4]] #M.Pst-Y.Pst
        
        #Select for the genes which are differentially expressed in M.Pst-M.Mock and M.Pst-Y.Pst and remove genes differentially expressed in Y.Pst-Y.Mock. Using the inverse of fisher's method
        data$padj =  -2*(log(1-temp[[4]]$padj)+log(1-temp[[2]]$padj)+log(temp[[1]]$padj))
        if (SUGs == T){
                data$padj =  -2*(log(1-temp[[4]]$padj)+log(temp[[2]]$padj)+log(1-temp[[1]]$padj))
        }
        data$padj =  pchisq(data$padj,2*3) #Uses the inverse of fisher's method to basically discount samples which aren't: 1. Differentiall expressed in mature plants 2. differentially expressed in M.Pst compared to Y.Pst 3. Are not differentially expressed in Y.Pst compared to Y.Mock
        
        data$de = "NO" #de = differentially expressed
        if (SUGs == F){
                data$de[data$log2FoldChange > lfc.cut & data$padj < adjp.cut] <- "UP"
                data$de[data$log2FoldChange < -lfc.cut & data$padj < adjp.cut] = "DOWN"    
        } else{
                data$de[data$log2FoldChange > lfc.cut & data$padj < adjp.cut] <- "DOWN"
                data$de[data$log2FoldChange < -lfc.cut & data$padj < adjp.cut] = "UP"  
        }

        mycolors <- c("blue", "red", "black")
        names(mycolors) <- c("DOWN", "UP", "NO")
        data$gene_symbol = getGeneName(rownames(data))
        data$gene_symbol[is.na(data$gene_symbol)] = ""
        data = data[is.na(data$log2FoldChange)==F & is.na(data$padj)==F,]
        
        data$delabel = NA
        n= 15
        data$delabel[order(data$log2FoldChange * (data$padj < adjp.cut) *(-log10(data$padj)), decreasing = T)][c(1:n,(nrow(data)-(n-1)):nrow(data))] = data$gene_symbol[order(data$log2FoldChange * (data$padj < adjp.cut) *(-log10(data$padj)), decreasing = T)][c(1:n,(nrow(data)-(n-1)):nrow(data))]
        data = data[order(data$log2FoldChange * (-log10(data$padj)), decreasing = T),]
        
        if (de < 0){
                data = data[order(data$log2FoldChange * (-log10(data$padj)), decreasing = F),]
                data = rownames(data[data$de == "DOWN", ])
                
        } else if(de > 0){
                data = rownames(data[data$de == "UP", ])
        }
        
        
        return(data)
}

#save.image("ARR-RNA-seqimage.RData")
load("ARR-RNA-seqimage.RData")

