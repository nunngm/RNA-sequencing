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
library(ggrepel)
library(NMF)
library(GenAnalysis)


hmcol = hcl_palettes(palette = "Berlin") #Setting the colour palatte
for_pca <- rlog(allData, blind=T)
rlogMat <- assay(for_pca) # just making a matrix of the counts that have been corrected for over-dispersion in a "blind" fashion

distsRL <- dist(t(rlogMat)) # Computes a distance matrix (Euclidian Distance)
mat <- as.matrix(distsRL)  # Make sure it is a matrix

rownames(mat) <- colnames(mat) <-   with(colData(allData), paste(age, infection, hpi, sep=" : "))

hc <- hclust(distsRL,method = "average")  # performs hierarchical clustering
par(mar=c(7,4,4,2)+0.1)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(20)# picking our colours
p=heatmap.2(mat, Rowv=as.dendrogram(hc),
            symm=TRUE, trace="none",
            col = rev(hmcol), margin=c(13,9))
ggsave("heatmap.pdf",plot = p)

##Begin intense heatmapping here

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


sub_24 = c(9,3,12,6) #samples collected at 24 hpi
sub_12 =c(8,2,11,5) #samples at 12 hpi
sub_0 = c(7,1,10,4) #Samples at 0 hpi
model_graph = c(7:9,1:3,10:12,4:6) #in order M mock 0-24 hpi, M Pst 0-24 hpi, young mock 0-24 hpi, young pst 0-24 hpi
full = c(7,1,10,4,8,2,11,5,9,3,12,6)
byage = c(7,10,1,4,8,11,2,5,9,12,3,6)
full = c(7:12,1:6)
pst = c(10:12,4:6)

hm_mat_DGEgenes = log_norm_counts[names(sig),sub_24] #Select only genes that were significantly differentially expressed at the for the indicated subgroup (ie 24 hpi)
rownames(hm_mat_DGEgenes) = objectSymbol[names(sig)] #use object symbol to replace accession #s with gene names where possible
hm_mat_DGEgenes =hm_mat_DGEgenes[rowSums(hm_mat_DGEgenes)!=0,]

hmcol = diverge_hcl(20,h = c(240,70),l = c(100,0),c=250,power = 3) #Selecting colours for the hierarchical cluster

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
hm_mat_DGEgenes = log_norm_counts[rownames(log_norm_counts) %in% gene,sub_24] #Uses dyplr to select rows of log_norm_counts where were found in the gene list created above
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
aheatmap(hm_mat_DGEgenes, Rowv = T, Colv = T, distfun= "euclidean",treeheight = 20,hclustfun = "average", scale = "row", cellwidth = 40)

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




