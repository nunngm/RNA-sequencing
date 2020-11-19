#Finally battling the time course hell

ageinf = factor(paste0(age,infection), levels = c("ymg", "ypst", "mmg", "mpst"))
infection = factor(infection, levels = c("mg","pst"))
age = factor(age, levels = c("y", "m"))
hpi = c(rep(0.25,3),rep(12,3),rep(24,3))
hpi = c(hpi,hpi,hpi,hpi)
hpi = as.numeric(hpi)
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
model_full = model.matrix(~0+age:infection+age:infection:hpi)
model_full

rawData <- DESeqDataSetFromHTSeqCount(design_full,design=model_full)

#Make group fator and trim any row that has less than 10 counts total
## This allows us to pick our treatment group (ie ymg0h rather than all of each type (all of 0 hpi or all of youbg or all of mock))

keep <- rowSums(counts(rawData)) >= 360
keep = rowMeans(counts(rawData)[, 1:18]) >=10 &rowMeans(counts(rawData)[, 19:36]) >=10
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

allData <- DESeq(rawData, test = "LRT", reduced =~ageinf+hpi)
resultsNames(allData)
res = results(allData, pAdjustMethod = "BH", alpha = 0.05)

#playing

res = results(allData,contrast = c(0,-0.5,-0.5,1,0,-5,0,1), pAdjustMethod = "BH", alpha = 0.05) #This works
#model #2
res = results(allData,contrast = c(0,1,0,-1,-1,-1,1), pAdjustMethod = "BH", alpha = 0.05)
#model #3
res = results(allData,contrast = c(0,1,-1,-1,1,0,1), pAdjustMethod = "BH", alpha = 0.05)
res = results(allData,contrast = c(0,0,-1,0,1,0,1), pAdjustMethod = "BH", alpha = 0.05)
#model #4
res = results(allData,contrast = c(0,-1,-1,1,-1,-1,-1,1), pAdjustMethod = "BH", alpha =0.05)
res = results(allData,contrast = c(0,-1,-1,1,0,-1,-1,1), pAdjustMethod = "BH", alpha =0.05)
