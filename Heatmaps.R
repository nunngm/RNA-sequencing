#Heat mapping
library(DESeq2)
library(RCurl)
library("RColorBrewer")
library(gplots)
library(topGO) 
library(Rgraphviz) 
library(dplyr)
library(ggplot2)
library(colorspace)
library(scales)
library(ggrepel)
library(NMF)
library(GenAnalysis)

#install packages
install.packages("devtools")
library(devtools)
install_github("nickytong/GenAnalysis")

# Let's reorganize this mess

## Colour palettes
### scale-white-salmon colour palette
hmcol1 = sequential_hcl(n = 10, h = 17, c = 117, l = c(64,100), power = 0.5, rev = T)
show_col(hmcol1)
hmcol2 = sequential_hcl(n = 10, h = 215, c = 41, l = c(46,100), power =0.5)
show_col(hmcol2)
hmcol = c(hmcol2[1:9], "#FFFFFF", hmcol1[2:10]) 
show_col(hmcol)

### Red-White-Blue colour palette
hmcol1 = sequential_hcl(n = 10, h = 13, c = 177, l = c(54,100), power = 0.4, rev = T)
show_col(hmcol1)
hmcol2 = sequential_hcl(n = 10, h = 266, c = 131, l = c(33,100), power =0.4)
show_col(hmcol2)
hmcol = c(hmcol2[1:9], "#FFFFFF", hmcol1[2:10]) 
show_col(hmcol)

#Bring in genes you want an expression heatmap for
mydata= read.table(file= "clipboard",sep= "\t",header =T)

selectGeneHeatmap = function(genesOfInterest,filename, colours =c("#0000FF","#3230EF","#4442E7", "#5251E3","#605FE1","#6E6DE1","#7E7DE3","#9191E7","#ABABED","#FFFFFF","#FFACAB","#FF918F","#FF7D7A","#FF6C68","#FF5C57","#FF4D47","#FF3E35","#FF2D20","#FF0000"),
                             padj = 0.05, width = 5, height =4, graph = F){
  comp = lapply(levels(hpi), function(x){ #Picks makes the comparisons for y.pst, m.mock, m.pst compared to y.mock at each time point
    list(results(allData, contrast = c("group", paste0("ypst", x), paste0("ymg", x)), alpha = padj, pAdjustMethod="BH", tidy = T), 
         results(allData ,contrast = c("group", paste0("mmg", x), paste0("ymg", x)), alpha = padj, pAdjustMethod="BH", tidy = T), 
         results(allData,contrast = c("group", paste0("mpst", x),paste0("ymg",x)), alpha = padj, pAdjustMethod="BH", tidy = T))
  })

  GOI.bool = comp[[1]][[1]]$row %in% genesOfInterest$accession #makes a list of T/F values that is just the genes of interest
  # make a matrices of log2-fold changes and adjusted p-values
  l2fc = cbind(comp[[1]][[1]]$log2FoldChange, comp[[1]][[2]]$log2FoldChange, comp[[1]][[3]]$log2FoldChange,
                   comp[[2]][[1]]$log2FoldChange, comp[[2]][[2]]$log2FoldChange, comp[[2]][[3]]$log2FoldChange,
                   comp[[3]][[1]]$log2FoldChange, comp[[3]][[2]]$log2FoldChange, comp[[3]][[3]]$log2FoldChange)
  pvals = cbind(comp[[1]][[1]]$padj, comp[[1]][[2]]$padj, comp[[1]][[3]]$padj,
                comp[[2]][[1]]$padj, comp[[2]][[2]]$padj, comp[[2]][[3]]$padj,
                comp[[3]][[1]]$padj, comp[[3]][[2]]$padj, comp[[3]][[3]]$padj)
  #Name the columns appropriately
  colnames(l2fc) = c("Y.Pst", "M.Mock", "M.Pst","Y.Pst", "M.Mock", "M.Pst", "Y.Pst", "M.Mock", "M.Pst" )
  colnames(pvals) = c("Y.Pst", "M.Mock", "M.Pst","Y.Pst", "M.Mock", "M.Pst", "Y.Pst", "M.Mock", "M.Pst")
  
  #subset to just genes of interest
  l2fc = l2fc[ GOI.bool, ]
  pvals = pvals[GOI.bool, ]
  
  #set the rownames appropriately (they have been reorganized by boolean selection)
  rownames(l2fc) = comp[[1]][[1]]$row[GOI.bool]
  rownames(pvals) = comp[[1]][[1]]$row[GOI.bool]
  
  #Rorder rows to be in order supplied in the genes of interest input
  l2fc = l2fc[unlist2(lapply(genesOfInterest$accession, function(x){grep(x, rownames(l2fc))})), ] #have to reorder
  
  #Set boxes with an adjusted p-val <0.05 to 0 (no sig diff)
  l2fc[pvals>0.05] = 0
  if(graph ==T){
    svg(filename = paste0(filename,".svg"),width = width, height = height)
    aheatmap(l2fc, color = colours, border_color = "#888888", Rowv = NA, Colv = NA, distfun= "euclidean",hclustfun = "average", scale = "none", labRow = genesOfInterest$label, breaks =-0.33)
    dev.off()
  }else{
    aheatmap(l2fc, color = pal(25), border_color = "#888888", Rowv = NA, Colv = NA, distfun= "euclidean",hclustfun = "average", scale = "none", labRow = genesOfInterest$label, breaks =-0.33)
  }
}



pal = colorRampPalette(c("blue", "white", "red"))
hmcol = hcl_palettes(palette = "Berlin") #Setting the colour palatte
for_pca <- rlog(allData, blind=F)
rlogMat <- assay(for_pca) # just making a matrix of the counts that have been corrected for over-dispersion in a "blind" fashion

distsRL <- dist(t(rlogMat)) # Computes a distance matrix (Euclidian Distance)
mat <- as.matrix(distsRL)  # Make sure it is a matrix


infection = c(rep("Mock",9),rep("Pst",9),rep("Mock",9),rep("Pst",9))
infection <- as.factor(infection)
length(infection)
rownames(mat) <- colnames(mat) <-   with(colData(allData), paste(age, infection, paste0(hpi, "h"), sep=":"))

hc <- hclust(distsRL
     ,method = "average"
     )  # performs hierarchical clustering
par(mar=c(7,4,4,2)+0.1)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)# picking our colours

tiff(filename = "heatmap.tiff", height = 1080, width = 1280) #opens a tiff device
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=T, trace="none",
          cexRow = 2,
          cexCol = 2,
          col = rev(hmcol), margin=c(13,9)) #prints to tiff
dev.off() #closes the file

##Begin intense heatmapping here
#Code for figure 1
allComp = compare.group()
thresh = 0.05
lfc = 1
genes_of_interest = unique(c(rownames(allComp$ypst0ymg0[allComp$ypst0ymg0$padj < thresh & abs(allComp$ypst0ymg0$log2FoldChange) > lfc, ]), rownames(allComp$ypst12ymg12[allComp$ypst12ymg12$padj < thresh & abs(allComp$ypst12ymg12$log2FoldChange) > lfc, ]), rownames(allComp$ypst24ymg24[allComp$ypst24ymg24$padj < thresh & abs(allComp$ypst24ymg24$log2FoldChange) > lfc, ]), rownames(allComp$mpst0mmg0[allComp$mpst0mmg0$padj < thresh & abs(allComp$mpst0mmg0$log2FoldChange) > lfc, ]), rownames(allComp$mpst12mmg12[allComp$mpst12mmg12$padj < thresh & abs(allComp$mpst12mmg12$log2FoldChange) > lfc, ]), rownames(allComp$mpst24mmg24[allComp$mpst24mmg24$padj < thresh & abs(allComp$mpst24mmg24$log2FoldChange) > lfc, ]) ))

genes_of_interest = unique(c(rownames(allComp$ypst0ymg0[allComp$ypst0ymg0$padj < thresh & abs(allComp$ypst0ymg0$log2FoldChange) > lfc, ]), rownames(allComp$ypst12ymg12[allComp$ypst12ymg12$padj < thresh & abs(allComp$ypst12ymg12$log2FoldChange) > lfc, ]), rownames(allComp$ypst24ymg24[allComp$ypst24ymg24$padj < thresh & abs(allComp$ypst24ymg24$log2FoldChange) > lfc, ]) ))

genes_of_interest = unique(c( rownames(allComp$mpst0mmg0[allComp$mpst0mmg0$padj < thresh & abs(allComp$mpst0mmg0$log2FoldChange) > lfc, ]), rownames(allComp$mpst12mmg12[allComp$mpst12mmg12$padj < thresh & abs(allComp$mpst12mmg12$log2FoldChange) > lfc, ]), rownames(allComp$mpst24mmg24[allComp$mpst24mmg24$padj < thresh & abs(allComp$mpst24mmg24$log2FoldChange) > lfc, ]) ))

mpstcomp = list(results(allData,contrast = c("group", "mpst0", "ymg0"), alpha = 0.05, pAdjustMethod="BH", tidy = T), results(allData,contrast = c("group", "mpst12", "ymg12"), alpha = 0.05, pAdjustMethod="BH", tidy = T), results(allData,contrast = c("group", "mpst24", "ymg24"), alpha = 0.05, pAdjustMethod="BH", tidy = T))

GOI = rownames(allComp$ypst0ymg0) %in% genes_of_interest
GOI = rownames(allComp$ypst0ymg0) %in% mydata$accsession

finalLFC = cbind(allComp$ypst0ymg0[GOI, 2], allComp$mmg0ymg0[GOI,2], mpstcomp[[1]][GOI, 3], allComp$ypst12ymg12[GOI, 2], allComp$mmg12ymg12[GOI, 2], mpstcomp[[2]][GOI, 3], allComp$ypst24ymg24[GOI, 2], allComp$mmg24ymg24[GOI, 2], mpstcomp[[3]][GOI, 3])

rownames(finalLFC) = rownames(allComp$ypst0ymg0)[rownames(allComp$ypst0ymg0) %in% genes_of_interest]
colnames(finalLFC) = c("Y.Pst", "M.Mock", "M.Pst","Y.Pst", "M.Mock", "M.Pst", "Y.Pst", "M.Mock", "M.Pst" )
finalLFC[ rownames(finalLFC)==names(labels),]

#Working on


#worked
hmcol1 = sequential_hcl(n = 10, h = 17, c = 117, l = c(64,100), power = 0.5, rev = T)
show_col(hmcol1)
hmcol2 = sequential_hcl(n = 10, h = 215, c = 41, l = c(46,100), power =0.5)
show_col(hmcol2)
hmcol = c(hmcol2[1:9], "#FFFFFF", hmcol1[2:10]) 
show_col(hmcol)

hmcol1 = sequential_hcl(n = 10, h = 13, c = 177, l = c(54,100), power = 0.4, rev = T)
show_col(hmcol1)
hmcol2 = sequential_hcl(n = 10, h = 266, c = 131, l = c(33,100), power =0.4)
show_col(hmcol2)
hmcol = c(hmcol2[1:9], "#FFFFFF", hmcol1[2:10]) 
show_col(hmcol)

heatmap_full <- ggplot(data = as.data.frame(finalLFC), 
                        aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill =value), colour="darkgrey") +
    scale_fill_gradient2(midpoint = 0, low = "#00798C", #mid = "white", 
                         high = "#FF6F59", na.value = "white", name = "Correlation", 
                         limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) + 
    labs(x="Ecotype and Condition", y="Clusters") +  theme_bw() + #theme(axis.text.y=element_blank(), plot.margin = margin(-0.75, 0, 0,0 , "cm")))
    theme(axis.text.y=element_text(angle=20,vjust=0.5), axis.text.x = element_text(angle = -30, vjust = 0.5, hjust = 0, size = 9) , plot.margin = margin(-0.75, 0, 0,0 , "cm"))
heatmap_full

mydata= read.table(file= "clipboard",sep= "\t",header =T)
labels = mydata$label
names(labels) = mydata$accsession

finalLFC = finalLFC[unlist2(lapply(names(labels), function(x){grep(x,rownames(finalLFC))})),]
rownames(finalLFC) = 

grep()
tiff("rplot2.tiff", width = 4000, height = 2000)
svg(filename = "ETI heatmap-reduv2.svg",width = 5, height = 4)
aheatmap(finalLFC, color = hmcol, border_color = "#888888", Rowv = NA, Colv = NA, distfun= "euclidean",hclustfun = "average", scale = "none", labRow = labels[rownames(finalLFC)], breaks =0)
dev.off()

pdf()
heatmap.2(finalLFC, Rowv=T, Colv = F,
          trace="none", col = hmcol,
          lhei = c(1,8), lwid = c(0.5,4),
          cexRow = 0.25, cexCol = 2,
          margin=c(10,9))
dev.off()



#Make a dataframe with the average normailized count number for each of the samples
allData = estimateSizeFactors(allData)
sizeFactors(allData) 
treatcounts = cbind(rowMeans(counts(allData,normalized = T)[,1:3]),
                    rowMeans(counts(allData,normalized = T)[,4:6]),
                    rowMeans(counts(allData,normalized = T)[,7:9]),
                    rowMeans(counts(allData,normalized = T)[,10:12]),
                    rowMeans(counts(allData,normalized = T)[,13:15]),
                    rowMeans(counts(allData,normalized = T)[,16:18]),
                    rowMeans(counts(allData,normalized = T)[,19:21]),
                    rowMeans(counts(allData,normalized = T)[,22:24]),
                    rowMeans(counts(allData,normalized = T)[,25:27]),
                    rowMeans(counts(allData,normalized = T)[,28:30]),
                    rowMeans(counts(allData,normalized = T)[,31:33]),
                    rowMeans(counts(allData,normalized = T)[,34:36]))
rownames(treatcounts) = rownames(counts(allData,normalized = T)) #Tags each row with appropriate accession number
colnames(treatcounts) = c("m:mg:0h","m:mg:12h","m:mg:24h","m:pst:0h","m:pst:12h","m:pst:24h","y:mg:0h","y:mg:12h","y:mg:24h","y:pst:0h","y:pst:12h","y:pst:24h")
treatcounts = as.data.frame(treatcounts)
write.csv(treatcounts,"meancounts.csv")

#Same thing as above but for the standard deviations for each of the treatment types
std = cbind(rowSds(counts(allData,normalized = T)[,1:3]),
            rowSds(counts(allData,normalized = T)[,4:6]),
            rowSds(counts(allData,normalized = T)[,7:9]),
            rowSds(counts(allData,normalized = T)[,10:12]),
            rowSds(counts(allData,normalized = T)[,13:15]),
            rowSds(counts(allData,normalized = T)[,16:18]),
            rowSds(counts(allData,normalized = T)[,19:21]),
            rowSds(counts(allData,normalized = T)[,22:24]),
            rowSds(counts(allData,normalized = T)[,25:27]),
            rowSds(counts(allData,normalized = T)[,28:30]),
            rowSds(counts(allData,normalized = T)[,31:33]),
            rowSds(counts(allData,normalized = T)[,34:36]))
std = as.data.frame(std)
colnames(std) = colnames(treatcounts)
rownames(std) = rownames(treatcounts)
write.csv(std,"stdcounts.csv")


##Log of the counts
log_norm_counts <- log2(treatcounts+ 1) #+1 to stop NAs from appearing incase a gene has 0 counts

#Takes a DESeq2 result data frame and turns it into a usable table
sigGenes = as.data.frame(res_shrunk_24h@listData)
rownames(sigGenes) = res_shrunk_24h@rownames

# sigGenes_12h = filter(sigGenes_12h, padj<0.05&abs(log2FoldChange)>0.58) #1.5 fold diff
sig = sigGenes[abs(sigGenes$log2FoldChange) > 1,]$padj #Genes that have a linear fold change >2 are kept
names(sig) <- rownames(sigGenes[abs(sigGenes$log2FoldChange) > 1,])
sig <- sig[complete.cases(sig)] #Removing incomplete rows (causes a ton of problems if kept)
sig = sig[sig<0.05] #Keep genes with p-val <0.05


sub_24 = c(9, 12, 3, 6) #samples collected at 24 hpi
sub_12 =c(8, 11, 2, 5) #samples at 12 hpi
sub_0 = c(7, 10, 1, 4) #Samples at 0 hpi
model_graph = c(7:9,1:3,10:12,4:6) #in order M mock 0-24 hpi, M Pst 0-24 hpi, young mock 0-24 hpi, young pst 0-24 hpi
full = c(7,1,10,4,8,2,11,5,9,3,12,6)
byage = c(7,10,1,4,8,11,2,5,9,12,3,6)
full = c(7:12,1:6)
pst = c(10:12,4:6)
##-----------------------------This is BROKEN and outputs uninterpreted results because the count # is logged and THEN scaled to the mean of the LOG counts this just screws with everything. To do this better it should use normalised counts and not scale the result.
#Additionally I would redo the analysis for genes so it was shrunk and pick all DEGs across all the time points then put them into 1 figure of Y.Mock, Y.Pst, M.Mock, M.Pst (Y.Pst-Y.Mock, M.Mock-Y.Mock, M.Pst-Y.Mock.)
hm_mat_DGEgenes = log_norm_counts[names(sig), sub_12] #Select only genes that were significantly differentially expressed at the for the indicated subgroup (ie 24 hpi)
rownames(hm_mat_DGEgenes) = unlist(lapply(names(sig), getGeneName))#use object symbol to replace accession #s with gene names where possible
hm_mat_DGEgenes =hm_mat_DGEgenes[rowSums(hm_mat_DGEgenes)!=0, ]
hm_mat_DGEgenes = hm_mat_DGEgenes[,2:4] - hm_mat_DGEgenes[,1]

tiff("rplot.tiff", width = 2000, height = 4000)
aheatmap(hm_mat_DGEgenes, color = hmcol, Rowv = T, Colv = T, distfun= "euclidean",treeheight = 20,hclustfun = "average", scale = "none", cellwidth = 40, breaks = 0)
     dev.off()

pdf()
heatmap.2(temp, Rowv=T, Colv = T,
     trace="none", col = hmcol,
     lhei = c(1,8), lwid = c(0.5,4),
     cexRow = 0.25, cexCol = 2,
           margin=c(10,9))
dev.off()

hmcol = diverge_hcl(21,h = c(240,70),l = c(100,0),c=250,power = 1.5) #Selecting colours for the hierarchical cluster
show_col(hmcol)

#Generating the heatmap
p=GenAnalysis::aheatmap(hm_mat_DGEgenes,color = hmcol, midpoint = 0,
                        clusterWithScaledData = T,scale="row", cluster_cols = F,
                        clustering_method = "average", clustering_distance_rows = "euclidean",
                        filename = "DEGs24hpi.pdf",border_color = NA , shrink = 1.5,returnTree = T, 
                        show_rownames = T, show_colnames = T,
                        cellwidth = 20, cellheight = 2,fontsize_row =   32
) 
#File output used as all parts of the image could be veiewed
#returnTree = T is set so that the genes in the same order (top to bottom) could be used to visually determine which genes were differentially epressed among the full collection
#This makes a heatmap using the colour scheme determined above
# This function was picked as it is available to cluster genes AFTER expression has been scaled thereby allowing to show genes which show a similar pattern of expression at very different levels of ecpression - > Scaling here works by determining the average expression across all treatment groups and scaling based on that


#A different heatmap which wasn't as effective
# NMF::aheatmap(hm_mat_DGEgenes, Rowv = T, Colv = NA,color = hmcol, distfun= "euclidean",treeheight = 100,hclustfun = "average", scale = "row", cellwidth = 10
#               #,cellheight = 5
#               ,filename = "gilly.pdf"
#               )


#Make a hierachical cluster based on a GO ID
gene = unique(gene_associations$DB_Object_ID[gene_associations$GO_ID=="GO:0009751"]) #A list of all genes associated with a specific GO ID
gene = gene[nchar(gene)==9] #Accesion #s are 9 character so anything not that length is booted

#Collect the log counts for genes of interest
hm_mat_DGEgenes = log_norm_counts[rownames(log_norm_counts) %in% gene,] #Uses dyplr to select rows of log_norm_counts where were found in the gene list created above
hm_mat_DGEgenes =hm_mat_DGEgenes[rowSums(hm_mat_DGEgenes)!=0,]
gene = rownames(hm_mat_DGEgenes)
rownames(hm_mat_DGEgenes) = objectSymbol[rownames(hm_mat_DGEgenes)]

#GenAnalysis Heatmap

#Set the color pallette for heat map
hmcol = diverge_hcl(20,h = c(240,70),l = c(90,0),c=250,power = 3)
p=GenAnalysis::aheatmap(hm_mat_DGEgenes,color = hmcol,midpoint = 0,cellwidth = 20,clusterWithScaledData = T,scale="row",cluster_cols = F,clustering_method = "average",clustering_distance_rows = "euclidean",filename = "silly.pdf",border_color = NA
                        , shrink = 1.5,returnTree = T,show_rownames = T,show_colnames = T #required for the next step
)

lab=gene[p$tree_row$order] #basically get the names of the genes and the order they are in the heatmap 

# NMF::aheatmap(hm_mat_DGEgenes,scale = "none", Rowv =T, Colv = NA, ##NMF heatmap
#          color = hmcol, #distfun= "euclidean",
#          treeheight = 80,#hclustfun = "average",
#          cellwidth = 20,
#        # ,filename="gilly.pdf" #,cellheight = 10
# )
#sigGenes = res_y

sigGenes = sigGenes[rownames(sigGenes)%in% lab,] #this is a handy way to collect  significant genes and rearrange them at the same time using dpylr
sigGenes$pvalue = as.integer( sigGenes$padj<0.05)
yet = lapply(lab, function(x){tmp = sigGenes[rownames(sigGenes)==x,]
return(tmp$pvalue)})
yet = unlist(yet)
yet[is.na(yet)] =0

hmcol = sequential_hcl(2,h =10,c=c(100,0),l = c(40,100),power = 1, rev =T)
cbind(1,yet)
GenAnalysis::aheatmap(cbind(1,yet),color = hmcol,cellwidth = 20,scale="none",cluster_rows = F,cluster_cols = F
                      # ,filename = "gillypValues.png",border_color = NA, show_rownames = F,show_colnames = F
                      , shrink = 1.5 #required for the next step
)

plotCounts(allData, 
           gene = "AT2G04550",
           intgroup="group",
           pch = 20, col = "red")

#NMF Heatmaps

#An old way of making a heatmap
log_norm_counts <- log2(counts(allData, normalized=TRUE) + 1)
res_mVa_12h <- results(allData, contrast = c("group", "mpst24", "mmg24"), alpha = 0.05, pAdjustMethod="BH")
summary( res_mVa_12h )
plotMA( res_mVa_12h , ylim =c(-5,5))
res_mVa_12h  <-  res_mVa_12h [order( res_mVa_12h$padj),]
gene = rownames( res_mVa_12h[1:20,])
DGEgenes = subset(res_mVa_12h,padj<0.05)
DGEgenes = rownames(subset(DGEgenes,abs( log2FoldChange)>2))
hm_mat_DGEgenes = log_norm_counts[DGEgenes,sub_24]

tiff("rplot.tiff", width = 2000, height = 4000)

# 3. Close the file
dev.off()
aheatmap(hm_mat_DGEgenes, Rowv = T, Colv = T, distfun= "euclidean",treeheight = 20,hclustfun = "average", scale = "none", cellwidth = 40)

p_fake <- rbeta(32833, 1,1) # you could also use runif(12627,1,1)
hist(p_fake, ylim=c (0,2500))
hist(res_mVa_12h$pvalue)

#calling the variables
goi= vector()
lfc = matrix(nrow =nrow(res),ncol = 6 )
agehpi = vector()
for (i in 1:nlevels(hpi)){ #I wanted to do y0 m0 y12 etc order so the 0 hours had to be grouped together
  for(j in nlevels(age):1){ #reversed the order so young then mature
    res = results(allData,contrast = c("group",paste0(levels(age)[j],"pst",levels(hpi)[i]),paste0(levels(age)[j],"mg",levels(hpi)[i])),alpha =0.05, pAdjustMethod = "BH")
    temp = as.data.frame(res@listData)
    rownames(temp) = res@rownames
    res = temp
    temp = res[abs(res$log2FoldChange)>1,]$padj #Selects genes which have a |l2fc|>1
    names(temp) = rownames(res[abs(res$log2FoldChange)>1,])
    temp = temp[complete.cases(temp)]
    temp = names(temp[temp<0.01]) #Selects genes with a <1% FDR
    goi=c(goi,temp)
    # names(temp) = rownames(res)
    lfc[,(i*2)-1+2-j] = res$log2FoldChange
    agehpi = c(agehpi,paste0(levels(age)[j],levels(hpi)[i]))
  }
}
rownames(lfc) = rownames(res)
colnames(lfc) = agehpi
goi = unique(goi)
goi = goi[order(goi)]
lfc = lfc[rownames(lfc)%in% goi,]
length(goi) ==nrow(lfc)
sum(rownames(lfc) ==goi) ==nrow(lfc)

hmcol = sequential_hcl(4,h = 65,c =150,l =c(100,0), power =2 ) #yellow
show_col(hmcol)
hmcol = c(sequential_hcl(4,h = 215,c =115,l =c(65,0), power =2 ) ,hmcol[(length(hmcol)-1):1])
show_col(hmcol)#Selecting colours for the hierarchical cluster
p=GenAnalysis::aheatmap(lfc,color = hmcol, midpoint = 0,
                        clusterWithScaledData = F,scale="none", cluster_cols = F,
                        clustering_method = "average", clustering_distance_rows = "euclidean",
                        filename = "bad.pdf",border_color = NA , shrink = 1.5,returnTree = T, 
                        show_rownames = F, show_colnames = T,
                        cellwidth = 20, truncate = T,Upper = 4, Lower = -4
)
rownames(lfc) = objectSymbol[rownames(lfc)]
write.table(lfc,file = "lfc.txt",sep = "\t")
#Used CLUSTER 3.0 pearson uncentered and average linkage
lfc2 = read.delim(file = "lfc2.txt",sep = "\t",stringsAsFactors = F)
lfc2 = lfc2[2:nrow(lfc2),1]
temp = lapply(rownames(lfc), function(x){tmp = grep(x,lfc2)
return(tmp)})
temp = unlist(temp)
length(temp) ==nrow(lfc)
lfc2 = lfc[temp,]
sum(rownames(lfc)[temp] == rownames(lfc2))
p=GenAnalysis::aheatmap(lfc2,color = hmcol, midpoint = 0,
                        clusterWithScaledData = F,scale="none", cluster_cols = F,
                        cluster_rows = F,
                        filename = "bad.pdf",border_color = NA , shrink = 1.5,returnTree = T, 
                        show_rownames = T, show_colnames = T,
                        cellwidth = 20,cellheight = 2,fontsize_row =   32, truncate = T,Upper = 4, Lower = -4
)




