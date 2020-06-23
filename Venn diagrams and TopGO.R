goi_y = y_up #sets the genes of interest for young plants
goi_m = m_up #sets the genes of interested for mature plants
goi_mm = mm_
goi_yy = yy_up

#Filter out genes that are not significant
sg_y = names(goi_y[goi_y<0.05])
sg_m = names(goi_m[goi_m<0.05])
sg_mm = names(goi_mm[goi_mm<0.05])
sg_yy =names(goi_yy[goi_yy<0.05])

#All genes sets all the unique genes in the environment and can be used to compare back too
allGenes = c(sg_y,sg_m)
allGenes = unique(allGenes)

#Set the colours for the venn diagram
colors = c("#c3db0f","#2cff21","#6b7fff","#ff4059","#de4dff")
venn.diagram(x = list(sg_y,sg_m,sg_mm),
             category.names = c("Y.Pst>Y.Mock","M.Pst>M.Mock","M.Mock>Y.Mock"),
             filename = "3venn12h.tiff",
             output = T,
             imagetype = "tiff",
             scaled = F,
             col = "black",
             fill = colors[1:3],
             cat.col = "black",
             cat.cex = 2,
             cat.dist = c(0.12, 0.12,0.04),
             margin =0.15)
#Two-way venn diagram - More useful
venn.diagram(x = list(sg_y,sg_m),
             category.names = c("Y.Mo>Y.Ps","M.Mo>M.Ps"),
             filename = "temp.tiff",
             output = T,
             inverted = T,
             rotation.degree = 180, #In this program the labels are stationary while the venn diagram is not. The circle with more genes will be put on the left therefore we have to rotate it so that the labels match the circle
             imagetype = "tiff",
             scaled = F,
             col = "black",
             fill = colors[1:2],
             cat.col = "black",
             cat.cex = 2,
             cat.dist = c(0.16, 0.16),
             margin =0.15)

#Specifies which part of the Venn diagram you would like to complete GO enrichment on
geneList = as.integer(allGenes %in% sg_y)
names(geneList) = allGenes
geneList[geneList==1] = as.integer( names(geneList[geneList==1])%in% sg_m)
geneList[geneList==1] = as.integer( names(geneList[geneList==1])%in% sg_mm)
sum(geneList)

geneList = as.factor(geneList) #This has to be a factor for TopGO

#Instead of just comparing to an environment of DEGs this allow comparison to all genes (which have an average expression >10 reads).
totalGenes = unique(c(names(goi_y),names(goi_m)))
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
showSigOfNodes(GOdata, score(fisher_GO), firstSigNodes = 15, useInfo = 'all')

#GO term -> significant genes in goi list
allGO = genesInTerm(GOdata)
sigGenes = lapply(allGO,function(x) x[x %in% names(geneList[geneList==1])] )
objectSymbol[sigGenes[["GO:0009725"]]]

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

venn = list(venn_y,venn_both,venn_m)
venn =data.frame(lapply(venn, "length<-", max(lengths(venn))))
venn = cbind(as.character(venn[,1]),as.character(venn[,1]),as.character(venn[,2]),as.character(venn[,2]),as.character(venn[,3]),as.character(venn[,3]))
venn = cbind(venn[,1:2],desvec[venn[,1]],venn[,3:4],desvec[venn[,3]],venn[,5:6],desvec[venn[,5]])
venn[,2] = objectSymbol[venn[,2]]
venn[,5] = objectSymbol[venn[,5]]
venn[,8] = objectSymbol[venn[,8]]
colnames(venn) = c("Young only","Gene name","Description","Both","Gene name","Description","Mature only","Gene name","Description")
write.csv(venn,"temp.csv")