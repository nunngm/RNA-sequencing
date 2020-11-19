#Set the
hp = "24"
res = compare.group(hpi= "24")
res[[2]]

y_up = res_y[res_y$log2FoldChange > 0,]$padj #collect genes where fold change is positive (up-regulated)
names(y_up) <- rownames(res_y[res_y$log2FoldChange > 0,]) #collect the genenames
y_up <- y_up[complete.cases(y_up)] #remove any NAs


goi_y = y_up #sets the genes of interest for young plants
goi_m = m_up #sets the genes of interested for mature plants
goi_mock = mock_up
goi_pst = pst_up

#Filter out genes that are not significant
sg_y = names(goi_y[goi_y < 0.01])
sg_m = names(goi_m[goi_m < 0.01])
sg_mock = names(goi_mock[goi_mock < 0.01])
sg_pst = names(goi_pst[goi_pst < 0.05])

#All genes sets all the unique genes in the environment and can be used to compare back too
allGenes = c(sg_y,sg_m
               #,sg_mock
              # ,sg_pst
     )
allGenes = unique(allGenes)

#Set the colours for the venn diagram
colors = c("#c3db0f","#2cff21","#6b7fff","#ff4059","#de4dff")
show_col(colors)

temp = c(tmp)
#Two-way venn diagram - More useful

vennGenes = list(microGenes,mature)
allGenes = unique(c(vennGenes[[1]], vennGenes[[2]]))
venn.diagram(x = list(vennGenes[[1]],vennGenes[[2]]),
             category.names = c("pen3 hyperinduced", "Mature uniquely up"),
             filename = "temp.tiff",
             output = T,
             inverted = T,
             rotation.degree = -180, #In this program the labels are stationary while the venn diagram is not. The circle with more genes will be put on the left therefore we have to rotate it so that the labels match the circle
             imagetype = "tiff",
             scaled = F,
             col = "black",
             fill = colors[1:2],
             cat.col = "black",
             cat.cex = 1,
             cat.dist = c(0.16, 0.16),
             margin =0.15)

#Specifies which part of the Venn diagram you would like to complete GO enrichment on
geneList = as.integer(allGenes %in% vennGenes[[1]])
names(geneList) = allGenes
geneList[geneList==1] = as.integer( names(geneList[geneList==1])%in% vennGenes[[2]])
geneList[geneList==1] = as.integer( names(geneList[geneList==1])%in% sg_mock)
sum(geneList)

getGeneName(names(geneList[geneList==1]))

geneList = as.factor(geneList) #This has to be a factor for TopGO

#Instead of just comparing to an environment of DEGs this allow comparison to all genes (which have an average expression >10 reads).
totalGenes = unique(rownames(res_y))
geneList2 = as.integer(totalGenes %in% names(geneList[geneList==1]))
names(geneList2) = totalGenes
sum(geneList2)
geneList2 = as.factor(geneList2)

#Does the GO enrichment (BP, MF, or CC)
GOdata <- new("topGOdata",
              ontology = "BP", 
              allGenes = geneList2,  
              annotationFun = annFUN.gene2GO, 
              gene2GO = gene_GO) 

fisher_test <- new("classicCount", testStatistic = GOFisherTest, name = "fisher_test")
fisher_GO <- getSigGroups(GOdata, fisher_test)
fisher_GO

#Lets see the most significant ones
table_GO <- GenTable(GOdata, Fisher = fisher_GO, topNodes = 80)
table_GO #A table of significant terms

#A tree of the significant terms
par(cex = 0.2)
showSigOfNodes(GOdata, score(fisher_GO), firstSigNodes = 30, useInfo = 'all')

#GO term -> significant genes in goi list
allGO = genesInTerm(GOdata)
sigGenes = lapply(allGO,function(x) x[x %in% names(geneList2[geneList2==1])] )
objectSymbol[sigGenes[["GO:0007165"]]]

tmp = data[rownames(data) %in% sigGenes[["GO:0007165"]], ]

for (i in 1:length(tmp)){
     df = tmp[[i]]
     for (j in 2:4){
          df[,j] = log2linear(df[,j])
     }
     if (i ==1){
          write.xlsx(df, "templax.xlsx", sheetName = as.character(i), col.names = T, row.names = F)
     } else {
          write.xlsx(df, "templax.xlsx", sheetName = as.character(i), col.names = T, row.names = F, append = T)
     }
}

df =  tmp[[3]][ tmp[[3]]$accession %in% temp, ]
for (i in 2:4){
     df[,i] = log2linear(df[,i])
}
write.xlsx(df, "vs.xlsx", sheetName = "RNASEQ", col.names = T, row.names = F)

df = final[unlist(lapply(final[,1],function(x){
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
})), ]
for ( i in 3:5){
     df[,i] = log2linear(df[,i])
}
write.xlsx( df, "vs.xlsx", sheetName = "Microarray",col.names = T, row.names = T, append = T
     )

#This outputs a file with all of the genes of interest from each part of a two-way venn diagram as a list with gene names and descriptions
venn_y = as.integer(allGenes %in% sg_y)
names(venn_y) = allGenes
venn_y[venn_y==1] = as.integer(! names(venn_y[venn_y==1])%in% sg_m)
venn_m = as.integer(allGenes %in% sg_m)
names(venn_m) = allGenes
venn_m[venn_m==1] = as.integer(! names(venn_m[venn_m==1])%in% sg_y)
venn_both = as.integer(allGenes %in% sg_y)
names(venn_both) = allGenes
venn_both[venn_both==1] = as.integer( names(venn_both[venn_both==1])%in% sg_m)
sum(venn_y)
sum(venn_both)
sum(venn_m)
venn_y = names(venn_y[venn_y==1])
venn_m = names(venn_m[venn_m==1])
venn_both = names(venn_both[venn_both==1])
# venn_y = cbind(venn_y, objectSymbol[venn_y],desvec[venn_y])
# venn = list(venn_y,venn_both,venn_m)

#prepping the excel file
df = genes.info(venn_y, hp)
# sum(is.na(df[,4]))
# df[,4] = log2linear(df[,4])
# df[,6] = log2linear(df[,6])
write.xlsx(df, "temp.xlsx", sheetName = "Y.Pst>Y.Mock 24 hpi", col.names = T, row.names = F)

df = genes.info(venn_both, hp)
write.xlsx(df, "temp.xlsx", sheetName = "Y.Pst>Y.Mock, M.Pst>M.Mock 24 hpi", col.names = T, row.names = F, append = T)

df = genes.info(venn_m, hp)
write.xlsx(df, "temp.xlsx", sheetName = "M.Pst>M.Mock 24 hpi", col.names = T, row.names = F, append = T)





column = c("Accession number", "Gene names", "Mean count", "Y.Pst-Y.Mock log2 fold-change", "q-val", "M.Pst-M.Mock log2 fold-change", "q-val", "Gene description")
df = cbind(venn_y, 
          objectSymbol[venn_y], 
          countMean[venn_y], 
          res_y$log2FoldChange[rownames(res_y) %in% venn_y], 
          res_y$padj[rownames(res_y) %in% venn_y], 
          res_m$log2FoldChange[rownames(res_m) %in% venn_y], 
          res_m$padj[rownames(res_m) %in% venn_y], 
          desvec[venn_y]
     )
df[is.na(df[,8]),8] = ""
colnames(df) = column
nrow(df)

write.xlsx(df, "temp.xlsx", sheetName = "Y.Pst>Y.Mock 12 hpi", col.names = T, row.names = F)

df = cbind(venn_both, 
          objectSymbol[venn_both], 
          countMean[venn_both], 
          res_y$log2FoldChange[rownames(res_y) %in% venn_both], 
          res_y$padj[rownames(res_y) %in% venn_both], 
          res_m$log2FoldChange[rownames(res_m) %in% venn_both], 
          res_m$padj[rownames(res_m) %in% venn_both], 
          desvec[venn_both]
     )
df[is.na(df[,8]),8] = ""
colnames(df) = column
nrow(df)

write.xlsx(df, "temp.xlsx", sheetName = "Y.Pst>Y.Mock, M.Pst>M.Mock 12 hpi", col.names = T, row.names = F, append = T)

df = cbind(venn_m, 
          objectSymbol[venn_m], 
          countMean[venn_m], 
          res_y$log2FoldChange[rownames(res_y) %in% venn_m], 
          res_y$padj[rownames(res_y) %in% venn_m], 
          res_m$log2FoldChange[rownames(res_m) %in% venn_m], 
          res_m$padj[rownames(res_m) %in% venn_m], 
          desvec[venn_m]
     )
df[is.na(df[,8]),8] = ""
colnames(df) = column
nrow(df)

write.xlsx(df, "temp.xlsx", sheetName = "M.Pst>M.Mock 12 hpi", col.names = T, row.names = F, append = T)

venn =data.frame(lapply(venn, "length<-", max(lengths(venn))))
venn = cbind(as.character(venn[,1]),as.character(venn[,1]),as.character(venn[,2]),as.character(venn[,2]),as.character(venn[,3]),as.character(venn[,3]))
venn = cbind(venn[,1:2],desvec[venn[,1]],venn[,3:4],desvec[venn[,3]],venn[,5:6],desvec[venn[,5]])
venn[,2] = objectSymbol[venn[,2]]
venn[,5] = objectSymbol[venn[,5]]
venn[,8] = objectSymbol[venn[,8]]
colnames(venn) = c("Young only","Gene name","Description","Both","Gene name","Description","Mature only","Gene name","Description")
write.csv(venn,"temp.csv")

#######--------------------------------------------------------------------------------------------
# Thre-way venn diagram code
colors = qualitative_hcl(3, palette = "warm",c = 90)
show_col(colors)
hour = c("0", "12", "24")
tmp = lapply(hour, find.unique)
temp = tmp
temp = list(tmp[[1]]$accession, tmp[[2]]$accession, tmp[[3]]$accession)
temp = list(microGenes, mature, young)

allGenes = c(temp[[1]], temp[[2]], temp[[3]])
allGenes = unique(allGenes)

venn.diagram(x = list(temp[[1]], temp[[2]], temp[[3]]),
             category.names = c("Microarray","M.Pst-M.Mock","Y.Pst-Y.Mock"),
             filename = "temp.tiff",
             output = T,
             imagetype = "tiff",
             euler.d = F,
             scaled = F,
             col = "black",
             fill = colors[1:3],
             cat.col = "black",
             cat.cex = 2,
             cat.dist = c(0.12, 0.12,0.04),
             margin =0.15)

venn = draw.triple.venn(427,144,3705,2,35,318,2,category = c("Y.Pst>Y.Mock","M.Pst>M.Mock","M.Mock>Y.Mock"),fill = colors,lty=1,cex=2,cat.cex=2,cat.col="black",cat.dist = c(0.12, 0.12,0.04),margin =0.15,euler.d = F,scaled = F)

#Specifies which part of the Venn diagram you would like to complete GO enrichment on
geneList = as.integer(allGenes %in% temp[[1]])
names(geneList) = allGenes
geneList[geneList==1] = as.integer(! names(geneList[geneList==1])%in% temp[[2]])
geneList[geneList==1] = as.integer( names(geneList[geneList==1])%in% temp[[3]])
sum(geneList)

getGeneName(names(geneList[geneList==1]))

geneList2 = as.factor(geneList)

#geneList = as.factor(geneList) #This has to be a factor for TopGO

#Instead of just comparing to an environment of DEGs this allow comparison to all genes (which have an average expression >10 reads).
totalGenes = unique(c(names(goi_y),names(goi_m),names(goi_mock)))
geneList2 = as.integer(totalGenes %in% names(geneList[geneList==1]))
names(geneList2) = totalGenes
sum(geneList2)
geneList2 = as.factor(geneList2)

#Does the GO enrichment (BP, MF, or CC)
GOdata <- new("topGOdata",
              ontology = "BP", 
              allGenes = geneList2,  
              annotationFun = annFUN.gene2GO, 
              gene2GO = gene_GO) 

fisher_test <- new("classicCount", testStatistic = GOFisherTest, name = "fisher_test")
fisher_GO <- getSigGroups(GOdata, fisher_test)
fisher_GO

#Lets see the most significant ones
table_GO <- GenTable(GOdata, Fisher = fisher_GO, topNodes = 50)
table_GO #A table of significant terms

#A tree of the significant terms
par(cex = 0.2)
showSigOfNodes(GOdata, score(fisher_GO), firstSigNodes = 20, useInfo = 'all')

#GO term -> significant genes in goi list
allGO = genesInTerm(GOdata)
sigGenes = lapply(allGO,function(x) x[x %in% names(geneList[geneList==1])] )
objectSymbol[sigGenes[["GO:0005634"]]]

objectSymbol[names(geneList[geneList==1])]

#_______________________ 4-way venn diagram___________________

colors = qualitative_hcl(4, palette = "warm",c = 90)
show_col(colors)
venn.diagram(x = list(sg_y,sg_m,sg_mock,sg_pst),
             category.names = c("Y.Pst>Y.Mock","M.Pst>M.Mock","M.Mock>Y.Mock","M.Pst>Y.Pst"),
             filename = "4way.tiff",
             output = T,
             imagetype = "tiff",
             euler.d = F,
             scaled = F,
             col = "black",
             fill = colors[1:4],
             cat.col = "black",
             cat.cex = 2,
             #cat.dist = c(0.12, 0.12,0.04),
             margin =0.15
     )

geneList = as.integer(allGenes %in% sg_m)
names(geneList) = allGenes
geneList[geneList==1] = as.integer(! names(geneList[geneList==1])%in% sg_y)
geneList[geneList==1] = as.integer( names(geneList[geneList==1])%in% sg_mock)
geneList[geneList==1] = as.integer( names(geneList[geneList==1])%in% sg_pst)
sum(geneList)
objectSymbol[names(geneList[geneList==1])]


venn.diagram(x = list(sg_y,sg_m,sg_mock),
             category.names = c("Y.Pst>Y.Mock","M.Pst>Y.Mock","M.Mock>Y.Mock"),
             filename = "temp.tiff",
             output = T,
             imagetype = "tiff",
             euler.d = F,
             scaled = F,
             col = "black",
             fill = colors[1:3],
             cat.col = "black",
             cat.cex = 2,
             cat.dist = c(0.12, 0.12,0.04),
             margin =0.15)