#Finally battling the time course hell

ageinf = factor(paste0(age,infection), levels = c("ymg", "ypst", "mmg", "mpst"))
infection = factor(infection, levels = c("mg","pst"))
age = factor(age, levels = c("y", "m"))
hpi = c(rep(1,3),rep(2,3),rep(3,3))
hpi = c(hpi,hpi,hpi,hpi)
hpi = as.integer(hpi)
length(hpi)
design_full <- data.frame(sample=names(files),
                          file=files,
                          age=age,
                          infection =infection,
                          hpi=hpi,
                          ageinf
)

design_full

##very general model that just has all the sample information
model_full <- formula(~ageinf+hpi+ageinf:hpi)
model_full = model.matrix(~0+ageinf + hpi + ageinf:hpi)#
model_full = model.matrix(~ageinf+ hpi + ageinf:hpi)# Model #1 this one mostly works although doesn't consider the diff between y.mock and m.mock
model_full = model.matrix(~age+infection+hpi+ageinf:hpi)#* Model #2
model_full = model.matrix(~age+infection:hpi+age:infection+age:infection:hpi)#* Model #3
model_full = model.matrix(~age+hpi+age:infection) #model 4
model_full = model.matrix(~0+age+infection+hpi+age:infection:hpi) 
model_full = model.matrix(~age+ infection + hpi + age:infection:hpi)#model 5
model_full = model.matrix(~0+age+infection:hpi+age:infection+age:infection:hpi)
model_full
model_full <- formula(~0+age:infection:hpi)#6
model_full = formula(~ageinf+hpi) #7

rawData <- DESeqDataSetFromHTSeqCount(design_full,design=model_full)

#Make group fator and trim any row that has less than 10 counts total
## This allows us to pick our treatment group (ie ymg0h rather than all of each type (all of 0 hpi or all of youbg or all of mock))

keep <- rowSums(counts(rawData)) >= 360
# keep = rowMeans(counts(rawData)[, 1:18]) >=10 &rowMeans(counts(rawData)[, 19:36]) >=10
rawData <- rawData[keep,]

#Lets redi the model design to encompass every sample individually
allData <- DESeq(rawData)

resultsNames(allData)

res = results(allData,contrast = c(0,-0.5,-0.5,1,0,-1,0,1), pAdjustMethod = "BH", alpha = 0.05) #This works!

res = res[order(res$log2FoldChange, decreasing = T), ]
res = res[is.na(res$padj)==F ,]
#res = res[ res$log2FoldChange > 0 & res$padj < 0.05,]
res2 = as.data.frame(res@listData)
rownames(res2) = as.character(res@rownames)
res2 = cbind(res2, getGeneName(rownames(res2)))
res2 = cbind(res2, objectSymbol[rownames(res2)])

# allData <- DESeq(rawData, test = "LRT", reduced =~ageinf+hpi)
# resultsNames(allData)
# res = results(allData, pAdjustMethod = "BH", alpha = 0.05)

#playing

res = results(allData,contrast = c(0,-0.5,-0.5,1,0,-5,0,1), pAdjustMethod = "BH", alpha = 0.05) #This works
#model #2
res = results(allData,contrast = c(0,1,0,-1,-1,-1,1), pAdjustMethod = "BH", alpha = 0.05)
#model #3
res = results(allData,contrast = c(0,1,-1,-1,1,0,1), pAdjustMethod = "BH", alpha = 0.05)
res = results(allData,contrast = c(0,0,-1,0,1,0,1), pAdjustMethod = "BH", alpha = 0.05)
res = results(allData,contrast = c(0,0,0,0,0,1,-1), pAdjustMethod = "BH", alpha = 0.05)
#model #4
res = results(allData, name = "agem", pAdjustMethod = "BH", alpha = 0.05)
res = results(allData,contrast = c(0,0,0,-1,1), pAdjustMethod = "BH", alpha = 0.05)
# model #5
res = results(allData,contrast = c(0,0,0,0,0,0,1,-1), pAdjustMethod = "BH", alpha = 0.05)
#6
res = results(allData,contrast = c(0,-1,0,1,0,-1,0,1,0,-1,0,1), pAdjustMethod = "BH", alpha = 0.05)
res = results(allData,contrast = c(-1,0,-1,1,-1,0,-1,1,-1,0,-1,1), pAdjustMethod = "BH", alpha = 0.05)
#7

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

model_full <- formula(~age+infection+hpi)
rawData <- DESeqDataSetFromHTSeqCount(design_full,design=model_full)
rawData$group <- factor(paste0(rawData$age,rawData$infection,rawData$hpi))
rawData$group = factor(rawData$group,levels=c("ymg0","ymg12","ymg24","ypst0","ypst12","ypst24","mmg0","mmg12","mmg24","mpst0","mpst12","mpst24"))
keep <- rowMeans(counts(rawData)) >= 10 #Genes which on average have less than 10 reads
#keep = rowMeans(counts(rawData)[, 1:18]) >=10 &rowMeans(counts(rawData)[, 19:36]) >=10
rawData <- rawData[keep,]

# Build DESeq object based on distinct treatment groups
rawData@design = ~0+group
allData <- DESeq(rawData)
resultsNames(allData)
res = results(allData, contrast = c(0,1,1, 0,-1,-1, 0,-1,-1, 0,1,1), pAdjustMethod = "BH", alpha = 0.05)
