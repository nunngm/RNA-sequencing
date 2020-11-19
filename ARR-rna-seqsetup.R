options(java.parameters = "-Xmx1024m")
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
library(scales)
library(xlsx)
library(stringi)

#laptop directory
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")

#Desktop directory
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")


gene_associations <- read.delim("gene_association_final.txt", comment.char = "!", header = FALSE, as.is = TRUE) 
colnames(gene_associations) <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
                                 "DB:Reference", "Evidence", "With_From", "Aspect", "DB_Object_Name",
                                 "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_by") 
# #I didn't want to look at locus ID's from TAIR so I took the first item of DB_Object_Synonym which is the gene ID and put it in the DB_Object_ID column as that is more useful
# gene_associations$DB_Object_ID = res
# gen_association_save <- sapply(gene_associations, function(x){paste(x, collapse = ", ")}) 
# gen_association_save <- data.frame("Gene Association" = gene_associations)
# write.table(gen_association_save, file = "gene_association_final.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = F)

# Trimming the dataframe so that it's only what we're interested in
gene_associations <- gene_associations[,c(2,3,5,7,9,14)]


##Make Gene-Go from scratch
# # Go through every unique gene and pull out any GO_ID associated with the given gene then give
# gene_GO <- lapply(unique(gene_associations$DB_Object_ID), function(x){tmp <- gene_associations %>% filter(DB_Object_ID == x)
# return(tmp$GO_ID)})
# names(gene_GO) <- unique(gene_associations$DB_Object_ID) 
# #Save it
# gene_GO_save <- sapply(gene_GO, function(x){paste(x, collapse = ", ")}) 
# gene_GO_save <- as.data.frame( gene_GO_save) # Making a dataframe
# write.table(gene_GO_save, file = "TAIR_to_GO.delim", sep = "\t", quote = FALSE, col.names = FALSE) 


#load Pre-meade Gene-to-go from file
gene_GO <- readMappings("TAIR_to_GO.delim")


files <- file.path("counts", list.files("counts"))


samples <- c("m_mg_0h_s1", "m_mg_0h_s2","m_mg_0h_s3","m_mg_12h_s1","m_mg_12h_s2","m_mg_12h_s3","m_mg_24h_s1","m_mg_24h_s2","m_mg_24h_s3","m_pst_0h_s1","m_pst_0h_s2","m_pst_0h_s3","m_pst_12h_s1","m_pst_12h_s2","m_pst_12h_s3","m_pst_24h_s1","m_pst_24h_s2","m_pst_24h_s3","y_mg_0h_s1","y_mg_0h_s2","y_mg_0h_s3","y_mg_12h_s1","y_mg_12h_s2","y_mg_12h_s3","y_mg_24h_s1","y_mg_24h_s2","y_mg_24h_s3","y_pst_0h_s1","y_pst_0h_s2","y_pst_0h_s3","y_pst_12h_s1","y_pst_12h_s2","y_pst_12h_s3","y_pst_24h_s1","y_pst_24h_s2","y_pst_24h_s3")
names(files) <- samples


sub_0 = c(1:3,10:12,19:21,28:30)
sub_12 = c(4:6,13:15,22:24,31:33)
sub_24 = c(7:9,16:18,25:27,34:36)

##Full data set

infection = c(rep("mg",9),rep("pst",9),rep("mg",9),rep("pst",9))
infection <- as.factor(infection)
length(infection)

hpi = c(rep(0,3),rep(12,3),rep(24,3))
hpi = c(hpi,hpi,hpi,hpi)
hpi = as.factor(hpi)
length(hpi)

age = c(rep("m",18),rep("y",18))
age = as.factor(age)
length(age)


design_full <- data.frame(sample=names(files),
                          file=files,
                          age=age,
                          infection =infection,
                          hpi=hpi
                )

design_full

##very general model that just has all the sample information
model_full <- formula(~age+infection+hpi)


rawData <- DESeqDataSetFromHTSeqCount(design_full,design=model_full)
rawData$infection = factor(rawData$infection, levels =  c("mg","pst"))
rawData$hpi = factor(rawData$hpi, levels = c("0","12","24"))
rawData$age = factor(rawData$age, levels = c("y","m"))

#Make group fator and trim any row that has less than 10 counts total
## This allows us to pick our treatment group (ie ymg0h rather than all of each type (all of 0 hpi or all of youbg or all of mock))
rawData$group <- factor(paste0(rawData$age,rawData$infection,rawData$hpi))
rawData$group = factor(rawData$group,levels=c("ymg0","ymg12","ymg24","ypst0","ypst12","ypst24","mmg0","mmg12","mmg24","mpst0","mpst12","mpst24"))

keep <- rowSums(counts(rawData)) >= 360
rawData <- rawData[keep,]

#Lets redi the model design to encompass every sample individually
rawData@design = ~group
allData <- DESeq(rawData)

##A really nice function (if a little messy) that completes the important comparisons (as I see it) by setting the constant 
#ie if you set the age to "y" then it will complete ypst - ymg at everytimepoint
# setting hpi restriects the function to pst - mg comparisons and m - y comparisons but wont compare mpst - ymg
compare.group = function(age = c("y", "m"), infection = c("mg", "pst"), hpi = c("0", "12", "24")){
     temp = lapply(age, function(x){tmp = paste0(x,infection)
                    return(tmp)
           })
     for (i in 1:length(temp)){
          temp[[i]] = lapply(temp[[i]], function(x){tmp = paste0(x, hpi)
                              return(tmp)
          })
     }
     out = character()
     for (i in 1:length(temp[[1]][[1]])){
          if(length(temp[[1]]) > 1){
               for (j in 1:length(temp)){
                    out = rbind(out, c(temp[[j]][[2]][i], temp[[j]][[1]][i]))
               }
          }
          if (length(temp) > 1){
               for (j in 1:length(temp[[1]])){
                    out = rbind(out, c(temp[[2]][[j]][i], temp[[1]][[j]][i]))
               }
          }
     }
     ##need to get rownames efficiently
     final = list()
     for (i in 1:nrow(out)){
          tmp = results(allData,contrast = c("group", out[i, 1], out[i, 2]), alpha = 0.05, pAdjustMethod="BH")
          temp = as.data.frame(tmp@listData)
          rownames(temp) = tmp@rownames
          remove(tmp)
          final[[i]] = temp
     }
     names(final) = paste0(out[,1], out[,2])
     return(final)
}

#Set-up the raw data comparison
#contrast = c(intgroup,Var1, Var2) -> contrasts Var 1 to reference (Var2)
#pAdjustMethod ="BH" -> is the false discovery rate method
## lfc Shrinkage doesn't matter here as p-values are calculated based on unshrunken values therefore shrinkage should only be used for data visualisation and for ranking of genes
res_y= results(allData,contrast = c("group", "ypst24", "ymg24"), alpha = 0.05, pAdjustMethod="BH")
temp = as.data.frame(res_y@listData)
rownames(temp) = res_y@rownames
res_y = temp

res_m = results(allData,contrast = c("group","mpst24","mmg24"),alpha =0.05, pAdjustMethod = "BH")
temp = as.data.frame(res_m@listData)
rownames(temp) = res_m@rownames
res_m = temp

res_mock = results(allData,contrast = c("group","mmg24","ymg24"),alpha =0.05, pAdjustMethod = "BH")
temp = as.data.frame(res_mock@listData)
rownames(temp) = res_mock@rownames
res_mock = temp

res_pst = results(allData, contrast = c("group", "mpst0", "ypst0"), alpha = 0.05, pAdjustMethod = "BH")
temp = as.data.frame(res_pst@listData)
rownames(temp) = res_pst@rownames
res_pst = temp

remove(temp)

##Collect up and down genes
# y_up = res_y[res_y$log2FoldChange > 0,]$padj #collect genes where fold change is positive (up-regulated)
# names(y_up) <- rownames(res_y[res_y$log2FoldChange > 0,]) #collect the genenames
# y_up <- y_up[complete.cases(y_up)] #remove any NAs
# 
# y_down = res_y[res_y$log2FoldChange < 0,]$padj #see above but for down-reg genes
# names(y_down) <- rownames(res_y[res_y$log2FoldChange < 0,]) 
# y_down <- y_down[complete.cases(y_down)]
# 
# m_up = res_m[res_m$log2FoldChange > 0,]$padj 
# names(m_up) <- rownames(res_m[res_m$log2FoldChange > 0,]) 
# m_up <- m_up[complete.cases(m_up)]
# 
# m_down = res_m[res_m$log2FoldChange < 0,]$padj 
# names(m_down) <- rownames(res_m[res_m$log2FoldChange < 0,]) 
# m_down <- m_down[complete.cases(m_down)]
# 
# mock_up = res_mock[res_mock$log2FoldChange > 0,]$padj #collect genes where fold change is positive (up-regulated)
# names(mock_up) <- rownames(res_mock[res_mock$log2FoldChange > 0,]) #collect the genenames
# mock_up <- mock_up[complete.cases(mock_up)]
# 
# mock_down = res_mock[res_mock$log2FoldChange > 0,]$padj #collect genes where fold change is positive (up-regulated)
# names(mock_down) <- rownames(res_mock[res_mock$log2FoldChange > 0,]) #collect the genenames
# mock_down <- mock_down[complete.cases(mock_down)]
# 
# pst_up = res_pst[res_pst$log2FoldChange > 0, ]$padj
# names(pst_up) = rownames(res_pst[res_pst$log2FoldChange > 0, ])
# pst_up = pst_up[complete.cases(pst_up)]
# 
# pst_down = res_pst[res_pst$log2FoldChange > 0, ]$padj
# names(pst_down) = rownames(res_pst[res_pst$log2FoldChange > 0, ])
# pst_down = pst_down[complete.cases(pst_down)]
# 


gene_descriptions = read.delim("gene_description.delim",sep = "\t",header = FALSE,stringsAsFactors = F,quote = "") 
gene_descriptions[,1]=tools::file_path_sans_ext(gene_descriptions[,1]) #removes ".1 or.2 ...etc" after each gene (these indicate the splice variant) 
gene_descriptions = gene_descriptions[,c(1,4)] #remove unnecessary columns
colnames(gene_descriptions) = c("accession", "description")
gene_descriptions = gene_descriptions %>% filter(!description == "") #remove splice variants with empty description
gene_descriptions = gene_descriptions %>% distinct(accession, .keep_all = T) #keep only the first splice variant with a description

#Makes a vector similar to object symbol to go from accession -> description
desvec = gene_descriptions[,2]
names(desvec) = gene_descriptions[,1] 

objectSymbol = lapply(unique(gene_associations$DB_Object_ID), function(x){tmp <- gene_associations %>% filter(DB_Object_ID == x)
return(tmp$DB_Object_Symbol)}) ##If gene assoc
names(objectSymbol) = unique(gene_associations$DB_Object_ID)
objectSymbol = unlist2(objectSymbol)
objectSymbol = objectSymbol[unique(names(objectSymbol))]

getGeneName = function(x){
     temp = objectSymbol[x]
     if (is.na(temp)==T){
          return(x)
     }
     return(temp)
}
countMean = res_y$baseMean
names(countMean) = rownames(res_y)

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

# view.gene = function(accession){
#   plotCounts(allData, 
#              gene = toupper(accession),
#              intgroup="group",
#              pch = 20, col = "red")
# }

# view.gene = function(accession){
#      geneCount <- plotCounts(allData, gene = toupper(accession), 
#                         intgroup = c("age", "infection","hpi"), returnData = TRUE)
#      geneCount$hpi <- as.numeric(as.character(geneCount$hpi))
#      geneCount$hpi[geneCount$hpi == 0] = 0.25
#      geneCount = cbind(geneCount, paste0(as.character(geneCount$age), as.character(geneCount$infection)))
#      colnames(geneCount)[ncol(geneCount)] = "ageinf"
#      
#      ggplot(geneCount,
#        aes(x = hpi, y = count, color = age, lty = infection, group = ageinf)) +
#           #geom_point() + #displays individual counts
#           stat_summary(fun=mean, geom="line", size = 1) + 
#           stat_summary(fun = mean,
#                     fun.min = function(x) {mean(x) - sd(x)}, 
#                     fun.max = function(x) {mean(x) + sd(x)}, 
#                     geom = "pointrange", lty =1 , size =1)+theme(panel.grid.major = element_blank(), 
#          panel.grid.minor = element_blank(),
#          panel.background = element_blank(), 
#          axis.line = element_line(colour = "black", size=1),
#          axis.title.x=element_text(size=15),
#          #axis.text.x=element_blank()),
#          axis.ticks=element_line(colour = "black", size =1),
#          axis.ticks.length = unit(5,"points") ,
#          axis.title.y = element_text(size=15),
#          legend.position = "right",
#          axis.text = element_text(size=15),
#          legend.text = element_text(size=15)
#      )
# }

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

     # ggsave(paste0(objectSymbol[toupper(accession)], ".png"), plot = p)
     # paste0(objectSymbol[toupper(accession)], ".png")
}

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


get.path = function(){
     return(paste(scan("clipboard", what = "string", sep = "\\"), collapse = "/"))
}

find.unique = function(hpi, reverse = F, down = F ){
     comp = compare.group(hpi = hpi)
     names(comp) = c("ypst_ymg", "mpst_mmg", "mmg_ymg", "mpst_ypst")
     temp = cbind( comp$ypst_ymg$log2FoldChange, comp$mmg_ymg$log2FoldChange, comp$mpst_mmg$log2FoldChange)
     pvals = cbind(comp$mmg_ymg$padj, comp$mpst_mmg$padj, comp$ypst_ymg$padj)
     
     # for (i in 1:nrow(pvals)){
     #      pvals[i, 1] = -2*sum(log(pvals[i, ]))
     #      pvals[i,1] = 1-pchisq(pvals[i,1], 2*length(pvals[i, ]))
     # }
     temp = cbind(rownames(comp$ypst_ymg), temp, pvals[ ,1], pvals[, 2], pvals[, 3])
     colnames(temp) = c("accession", "ypst_ymg", "mmg_ymg", "mpst_mmg", "qval", "qmpst", "qypst")
     temp = as.data.frame(temp)
     for (i in 2:ncol(temp)){
          temp[, i] = as.numeric(as.character(temp[, i]))
     }
     
     if (reverse == F & down == F){
          # temp = subset(temp, qval < 0.01 | qmpst < 0.01)
          temp = subset(temp, qmpst <0.1)
          temp = subset(temp, mpst_mmg > 0)
          temp = subset(temp, ypst_ymg < 0 | qypst > 0.1)
          temp = subset(temp, mmg_ymg + mpst_mmg - ypst_ymg > log2(1.5))
          temp = temp[order(temp$mmg_ymg + temp$mpst_mmg - temp$ypst_ymg, decreasing = T), ]
          temp = cbind(temp, objectSymbol[temp$accession], desvec[temp$accession])
     } else if (reverse == T & down == F){
          temp = subset(temp, ypst_ymg > 0 & qypst < 0.05)
          temp = subset(temp, mpst_mmg < 0 | qmpst > 0.05)
          temp = subset(temp, ypst_ymg - mmg_ymg - mpst_mmg > log2(1.5))
          temp = temp[order(temp$ypst_ymg - temp$mmg_ymg - temp$mpst_mmg, decreasing = T), ]
          temp = cbind(temp, objectSymbol[temp$accession], desvec[temp$accession])
     } else {
             temp = subset(temp, qmpst < 0.05)
             temp = subset(temp, mpst_mmg < 0)
             temp = subset(temp, ypst_ymg > 0 | qypst >0.05)
             temp = temp[order(temp$mmg_ymg + temp$mpst_mmg - temp$ypst_ymg, decreasing = F), ]
             temp = cbind(temp, objectSymbol[temp$accession], desvec[temp$accession])
     }

     return(temp)
}

find.volcano = function(hour, n = 15, de = 0) {
        temp = compare.group(hpi= hour)
        data = temp[[4]] #M.Pst-Y.Pst
        
        #Select for the genes which are differentially expressed in M.Pst-M.Mock and M.Pst-Y.Pst and remove genes differentially expressed in Y.Pst-Y.Mock. Using the inverse of fisher's method
        data$padj =  -2*(log(1-temp[[4]]$padj)+log(1-temp[[2]]$padj)+log(temp[[1]]$padj))
        data$padj =  pchisq(data$padj,2*3) #Uses the inverse of fisher's method to basically discount samples which aren't: 1. Differentiall expressed in mature plants 2. differentially expressed in M.Pst compared to Y.Pst 3. Are not differentially expressed in Y.Pst compared to Y.Mock
        
        data$de = "NO" #de = differentially expressed
        data$de[data$log2FoldChange > 0.585 & data$padj < 0.05] <- "UP"
        data$de[data$log2FoldChange < -0.585 & data$padj < 0.05] = "DOWN"
        mycolors <- c("blue", "red", "black")
        names(mycolors) <- c("DOWN", "UP", "NO")
        data$gene_symbol = getGeneName(rownames(data))
        data$gene_symbol[is.na(data$gene_symbol)] = ""
        data = data[is.na(data$log2FoldChange)==F & is.na(data$padj)==F,]
        
        data$delabel = NA
        n= 15
        data$delabel[order(data$log2FoldChange * (data$padj <0.05) *(-log10(data$padj)), decreasing = T)][c(1:n,(nrow(data)-(n-1)):nrow(data))] = data$gene_symbol[order(data$log2FoldChange * (data$padj <0.05) *(-log10(data$padj)), decreasing = T)][c(1:n,(nrow(data)-(n-1)):nrow(data))]
        data = data[order(data$log2FoldChange * (-log10(data$padj)), decreasing = T),]
        
        if (de < 0){
                data = data[order(data$log2FoldChange * (-log10(data$padj)), decreasing = F),]
                data = rownames(data[data$de == "DOWN", ])
                
        } else if(de > 0){
                data = rownames(data[data$de == "UP", ])
        }
        
        
        return(data)
}

# #Playing around with contrasts for timepoint analysis
# targets = as.data.frame(cbind(c("f1", "f2", "f3", "f4", "f5", "f6"), c("wt.0h", "wt.6h", "wt.24h", "mut.0h", "mut.6h", "mut.24h")))
# colnames(targets) = c("FileName", "Targets")
# lev <- c("wt.0h","wt.6h","wt.24h","mut.0h","mut.6h","mut.24h")
# f <- factor(targets$Target, levels=lev)
# design <- model.matrix(~0+f)
# colnames(design) <- lev
# fit <- lmFit(eset, design)
# 
# 
# cont.dif <- makeContrasts(
#          Dif6hr =(mu.6hr-mu.0hr)-(wt.6hr-wt.0hr),
#          Dif24hr=(mu.24hr-mu.6hr)-(wt.24hr-wt.6hr),
#          levels=design)
# fit2 <- contrasts.fit(fit, cont.dif)
# fit2 <- eBayes(fit2)
# topTable(fit2, adjust="BH")
# 
# 
# 
# files <- file.path("counts", list.files("counts"))
# samples <- c("m_mg_0h_s1", "m_mg_0h_s2","m_mg_0h_s3","m_mg_12h_s1","m_mg_12h_s2","m_mg_12h_s3","m_mg_24h_s1","m_mg_24h_s2","m_mg_24h_s3","m_pst_0h_s1","m_pst_0h_s2","m_pst_0h_s3","m_pst_12h_s1","m_pst_12h_s2","m_pst_12h_s3","m_pst_24h_s1","m_pst_24h_s2","m_pst_24h_s3","y_mg_0h_s1","y_mg_0h_s2","y_mg_0h_s3","y_mg_12h_s1","y_mg_12h_s2","y_mg_12h_s3","y_mg_24h_s1","y_mg_24h_s2","y_mg_24h_s3","y_pst_0h_s1","y_pst_0h_s2","y_pst_0h_s3","y_pst_12h_s1","y_pst_12h_s2","y_pst_12h_s3","y_pst_24h_s1","y_pst_24h_s2","y_pst_24h_s3")
# names(files) <- samples
# 
# ##Full data set
# infection = c(rep("mg",9),rep("pst",9),rep("mg",9),rep("pst",9))
# infection <- as.factor(infection)
# length(infection)
# 
# hpi = c(rep(0,3),rep(12,3),rep(24,3))
# hpi = c(hpi,hpi,hpi,hpi)
# hpi = as.factor(hpi)
# length(hpi)
# 
# age = c(rep("m",18),rep("y",18))
# age = as.factor(age)
# length(age)
# 
# 
# design_full <- data.frame(file=files,
#                           age=age,
#                           infection =infection,
#                           hpi=hpi,
#                           factor(paste(design_full$age, design_full$infection, design_full$hpi, sep = "_"))
# )
# temp = cbind(design_full, pas)
# colnames(design_full)[length(colnames(design_full))] = "treatment"
# design_full
# 
# 
# cont.totdif <- makeContrasts(
#         Dif12hr = ((m_pst_12-m_pst_0)-(m_mg_12-m_mg_0))-((y_pst_12-y_pst_0)-(y_mg_12-y_mg_0)),
#         Dif24hr = ((m_pst_24-m_pst_12)-(m_mg_24-m_mg_12))-((y_pst_24-y_pst_12)-(y_mg_24-y_mg_12)),
#         levels=design)
# fit2 <- contrasts.fit(fit, cont.dif)
# fit2 <- eBayes(fit2)
# topTable(fit2, adjust="BH")
