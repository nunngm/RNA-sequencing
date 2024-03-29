## Base file forked from caitsimop/eutrema_drought_transcriptomics/drought_wgcna.R
## want to extract log2fold changes for all these genes in each sample
## want to change the log2fold value to 0 if not significant!!
## Let's try a co-expression network first
library(WGCNA)
library(reshape2)
library(ggdendro)
library(ggplot2)
library(tidyverse)
library(viridis)
library(patchwork) #now on CRAN! huzzah
library(phylobase)
library(gplots)
library(RColorBrewer)
library(topGO)
library(WriteXLS)
library(zoo)
# # somethings depend on: ARR-rna-seqsetup.R
# # uncomment if you want to source this file
# setwd("C:/Users/Garrett/Documents/Local Git/RNA-sequencing")
# source("ARR-rna-seqsetup.R")
allgenes = counts(allData, normalized = T)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

gsg <- goodSamplesGenes(t(allgenes), verbose = 3);
gsg$allOK
# allgenes <- t(allgenes)
sampleTree = hclust(dist(t(allgenes)), method = "average")

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0));
par(mfrow = c(1,1))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
    cex.axis = 1.5, cex.main = 2)
dev.off()
## Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(allgenes, powerVector = powers, verbose = 5, 
                        networkType = "signed hybrid")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## 6 looks like a good fit

# need to do it block wise because my computer doesnt have enough RAM for this
network = blockwiseModules(t(allgenes), power = 6,
                           TOMType = "none",
                           networkType = "signed hybrid",
                           randomSeed = 31138, corType = "pearson", 
                           minModuleSize = 30,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "droughtTOM",
                           verbose = 3, maxBlockSize = 10000)

## if you save this, it is much quicker to run
save.image("limma.blockwiseNetwork.RData")
load("blockwiseNetwork.RData")


clustcolours <- labels2colors(network$colors)
MEs = moduleEigengenes(t(allgenes), clustcolours)$eigengenes




nGenes = nrow(allgenes)
nSamples = ncol(allgenes);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(t(allgenes), clustcolours)$eigengenes
MEs = orderMEs(MEs0)
ecotype_condition <- c(rep("M_Mock_0.25", 3), 
                       rep("M_Mock_12", 3),
                       rep("M_Mock_24", 3),
                       rep("M_Pst_0.25", 3),
                       rep("M_Pst_12", 3), 
                       rep("M_Pst_24", 3),
                       rep("Y_Mock_0.25", 3),
                       rep("Y_Mock_12", 3),
                       rep("Y_Mock_24", 3),
                       rep("Y_Pst_0.25", 3),
                       rep("Y_Pst_12", 3),
                       rep("Y_Pst_24", 3)
                       )

datTraits <- data.frame(ecocond = as.factor(ecotype_condition))
rownames(datTraits) <- colnames(allgenes)
traits <- model.matrix(~ ., data=datTraits,
                       contrasts.arg = lapply(datTraits, contrasts, contrasts=FALSE))
traits <- traits[,-1]
moduleTraitCor <- WGCNA::cor(MEs, traits, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,31);

pal <- colorRampPalette(c("#E7B800","#E7B800", "white", "#67A7C1", "#67A7C1"),
                        alpha=TRUE, interp="spline")(100)
droughtcol <- c("#C0C8C8","#8D9B9B",  "#6A90D1",
                "#C0C8C8","#8D9B9B",  "#6A90D1",
                "#C0C8C8","#8D9B9B",  "#6A90D1",
                "#C0C8C8","#8D9B9B",  "#6A90D1")

colnames(moduleTraitCor) <- c("M_Mock_0.25", "M_Mock_12", "M_Mock_24", "M_Pst_0.25", "M_Pst_12", "M_Pst_24", "Y_Mock_0.25", "Y_Mock_12", "Y_Mock_24", "Y_Pst_0.25", "Y_Pst_12", "Y_Pst_24")
sample_clust <- hclust(dist(t(moduleTraitCor), method="euclidean"))
module_clust <- hclust(as.dist(1-cor(t(moduleTraitCor), method="pearson")))

rownames(moduleTraitCor) <- rownames(moduleTraitCor) %>% substr(.,3,nchar(rownames(moduleTraitCor)))

head(moduleTraitCor)

colnames(moduleTraitCor) <- c("M_Mock_0.25", "M_Mock_12", "M_Mock_24", "M_Pst_0.25", "M_Pst_12", "M_Pst_24", "Y_Mock_0.25", "Y_Mock_12", "Y_Mock_24", "Y_Pst_0.25", "Y_Pst_12", "Y_Pst_24")

groupedByHpi = c(7,10,1,4,8,11,2,5,9,12,3,6)

# set non-significant correlations to 0
moduleTraitCor = moduleTraitCor * (moduleTraitPvalue < 0.05)

corr_clust <- moduleTraitCor[module_clust$order, groupedByHpi]
corr_clust <- moduleTraitCor[module_clust$order, sample_clust$order]
correlation_ggplot <- corr_clust%>% melt()
dx <- dendro_data(sample_clust)
dy <- dendro_data(module_clust)

# helper function for creating dendograms
ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

# x/y dendograms
px <- ggdend(dx$segments) 
py <- ggdend(dy$segments) + coord_flip()


clust.of.interest <- c("green", "greenyellow", "skyblue3", "saddlebrown", "grey60", "darkolivegreen", "navajowhite2", "salmon4", "magenta", "yellow")
clust.of

correlation_ggplot$Var2 = factor(correlation_ggplot$Var2, levels = c("Y_Mock_0.25", "Y_Pst_0.25", "Y_Mock_12", "Y_Pst_12", "Y_Mock_24", "Y_Pst_24", "M_Mock_0.25", "M_Pst_0.25", "M_Mock_12", "M_Pst_12", "M_Mock_24", "M_Pst_24"))
par(cex=0.5)
midpoint = 0, low = "#00798C", mid = "white", 
high = "#FF6F59"
(heatmap <- ggplot(data = correlation_ggplot %>% filter(., Var1 %in% clust.of.interest), 
                   aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill =value), colour="darkgrey") +
    scale_fill_gradient2(midpoint = 0, low = "#0000FF", mid = "white", 
                         high = "#FF0000", na.value = "white", name = "Correlation", 
                         limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) + 
    labs(x="", y="") +  theme_bw() + #theme(axis.text.y=element_blank(), plot.margin = margin(-0.75, 0, 0,0 , "cm")))
    theme(axis.text.y=element_text(angle=0,vjust=0.5, size = 15), axis.text.x = element_text(angle = -30, vjust = 1, hjust = 0, size = 9), plot.margin = margin(-0.75, 0, 0,0 , "cm")))
(dendro_heatmap <- px + heatmap + plot_layout(ncol = 1, heights = c(1, 5))) # + plot_layout(ncol=2, widths=(1,5))
(heatmap = heatmap + plot_layout(ncol = 1, heights = 1))

svg(filename = "select_clust_heatmap.svg", width = 8, height = 6)
dev.off()

ggsave("select_clust_heatmap.pdf", dendro_heatmap,
       height=6, width = 7, units = "in")

pdf(file = "all_dendo_clust_heatmap-limma.pdf")
heatmap_full <- ggplot(data = correlation_ggplot, 
                        aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill =value), colour="darkgrey") +
    scale_fill_gradient2(midpoint = 0, low = "#00798C", #mid = "white", 
                         high = "#FF6F59", na.value = "white", name = "Correlation", 
                         limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) + 
    labs(x="Ecotype and Condition", y="Clusters") +  theme_bw() + #theme(axis.text.y=element_blank(), plot.margin = margin(-0.75, 0, 0,0 , "cm")))
    theme(axis.text.y=element_text(angle=20,vjust=0.5), axis.text.x = element_text(angle = -30, vjust = 0.5, hjust = 0, size = 9) , plot.margin = margin(-0.75, 0, 0,0 , "cm"))
dendro_heatmap_full <- px + heatmap_full + plot_layout(ncol = 1, heights = c(1, 5))
dendro_heatmap_full = heatmap_full + plot_layout(ncol = 1, heights =  5)
dendro_heatmap_full
dev.off()

ggsave("./plots/all_clust_heatmap.pdf", dendro_heatmap_full,
       height=10, width = 8, units = "in")

geneClusters <- network$colors
genenames <- rownames(allgenes)
names(clustcolours) <- genenames

genesOI = names(clustcolours[clustcolours == "greenyellow"])
     objectSymbol[genesOI]
temp = cbind(genesOI, objectSymbol[genesOI])
View(temp)

volGenes = find.volcano(hour = "12")
sum(genesOI %in% rownames(volGenes[volGenes$de == "UP", ]))
selectGenes = genesOI[genesOI %in% rownames(volGenes[volGenes$de == "UP", ])]
selectGenes = objectSymbol[selectGenes]
selectGenes[selectGenes == "AT2G13810"]
# 308/848
name = "green"

genesOI = names(clustcolours[clustcolours == name])
temp = genes.info(genesOI, hpi = "0", mock = T)
write.xlsx(temp, paste0(name, ".xlsx"), sheetName = "0.25h", col.names = T, row.names = T)
temp = genes.info(genesOI, hpi = "12", mock = T)
write.xlsx(temp, paste0(name, ".xlsx"), sheetName = "12h", col.names = T, row.names = T, append = T)
temp = genes.info(genesOI, hpi = "24", mock = T)
write.xlsx(temp, paste0(name, ".xlsx"), sheetName = "24h", col.names = T, row.names = T, append = T)

## What about ALL DEGS?? 
## let's put all DEGs into a list, then are able to see which cluster(s) entire groups of degs are in!
## this depends on output from deseq2_drought.R

AUGs = unique(c(find.unique("0")$accession, find.unique("12")$accession, find.unique("24")$accession))
colOI = clustcolours[names(clustcolours) %>% toupper() %in% AUGs]

SUGs = unique(c(find.unique("0", reverse = T)$accession, find.unique("12", reverse = T)$accession, find.unique("24", reverse = T)$accession))
colOI = clustcolours[names(clustcolours) %>% toupper() %in% SUGs]

colOI = clustcolours[names(clustcolours) %>% toupper() %in% rownames(temp)]

data = compare.group(altH = "greater")
genesOI = rowMeans(cbind((-data[[1]]$stat+data[[2]]$stat-data[[3]]$stat+data[[4]]$stat)/sqrt(4), (data[[8]]$stat-data[[7]]$stat+data[[6]]$stat-data[[5]]$stat)/sqrt(4), (data[[12]]$stat-data[[11]]$stat+data[[10]]$stat-data[[9]]$stat)/sqrt(4)))
genesOI = genesOI >1.65
AUGs = rownames(data[[1]][ genesOI, ])

volc = list(find.volcano("0"), find.volcano("12"), find.volcano("24"))
uod = "UP"

volc = unique(c(rownames(volc[[1]][volc[[1]]$de == uod, ]), rownames(volc[[2]][volc[[2]]$de == uod, ]), rownames(volc[[3]][volc[[3]]$de == uod, ]) ))
colOI = clustcolours[names(clustcolours) %>% toupper() %in% volc]

# summarize findings from modules
tbl = as.data.frame(table(colOI))
tbl = tbl[order(tbl$Freq, decreasing = T), ]
tbl = cbind(tbl, t(as.data.frame(lapply(tbl$colOI, function(x) {tmp = length(clustcolours[clustcolours %in% x])}))))
tbl = cbind(tbl, c(tbl$Freq/tbl[, 3]))
colnames(tbl) = c("colOI", "freq", "size", "ratio")
tbl[order(tbl$ratio, decreasing = T), ]
tbl[order(tbl$freq, decreasing = T), ]

clustcolours[names(clustcolours) %>% toupper() %in% pattern4genes]


genesOI = clustcolours[clustcolours %in% "navajowhite2"]
objectSymbol[names(genesOI)]

##################
## GO enrichment #
##################

geneID2GO <- readMappings(file="./data_files/goannotation.txt", sep="\t")

GOterms <- read.csv("./data_files/goannotation.txt", sep="\t", header=FALSE)
colnames(GOterms) <- c("GeneID", "GOterms")
allGenes <- rownames(allgenes)


fisherStat <- new("classicCount", testStatistic = GOFisherTest, name =
                    "Fisher test")

### make gene lists for the clusters we're interested in:
clustersDEGS <- c("lightcyan1", "lightyellow", "purple", "pink", "blue", "darkslateblue", 
                  "turquoise", "coral1", "red")

genelist <- function(colour){
  genes <- names(clustcolours)[clustcolours == colour]
  genestosend <- fpkm[row.names(unique(fpkm[,c(1,10:52)])),]
  x <- genestosend[genestosend$gene_id %in% genes,] 
  return(x)
}


cluster_genes <- list()

for (x in clustersDEGS){
  cluster_genes[[x]] <- genelist(x)
}

#WriteXLS(cluster_genes,
#         ExcelFileName="DEGcluster_geneswithred.xls",
#         FreezeCol = 1, FreezeRow = 1)



goenrichment <- function(colour, allclusters, genenames, geneGOIDs){
  fisherStat <- new("classicCount", testStatistic = GOFisherTest, name =
                      "Fisher test")
  genes <- names(allclusters)[allclusters == colour]
  genes <- factor(as.integer(genenames %in% genes))
  names(genes) <- genenames
  GOdata <- new("topGOdata", ontology = "BP", allGenes = genes,
                annot = annFUN.gene2GO, gene2GO = geneGOIDs)
  resFisher <- getSigGroups(GOdata, fisherStat)
  pvals <- score(resFisher)
  adjustPval <- p.adjust(pvals, method = "fdr", n=length(pvals))
  sigGenes <- adjustPval[which(adjustPval <= 0.05)]
  results <- GenTable(GOdata, classic = resFisher,
                      topNodes = length(sigGenes))
  res <- merge(results, sigGenes, by.x="GO.ID", by.y=0)
  colnames(res)[7] <- "p.adj" 
  return(res)
}


brown_res <- goenrichment(colour="brown", allclusters = clustcolours, genenames = allGenes,
                          geneGOIDs = geneID2GO)


turquoise_res <- goenrichment(colour="turquoise", allclusters = clustcolours, genenames = allGenes,
                              geneGOIDs = geneID2GO)

blue_res <- goenrichment(colour="blue", allclusters = clustcolours, genenames = allGenes,
                         geneGOIDs = geneID2GO)


brown_res <-  goenrichment(colour="brown", allclusters = clustcolours, genenames = allGenes,
                           geneGOIDs = geneID2GO)

red_res <- goenrichment(colour="red", allclusters = clustcolours, genenames = allGenes,
                        geneGOIDs = geneID2GO)
yellow_res <- goenrichment(colour="yellow", allclusters = clustcolours, genenames = allGenes,
                           geneGOIDs = geneID2GO)
pink_res <- goenrichment(colour="pink", allclusters = clustcolours, genenames = allGenes,
                         geneGOIDs = geneID2GO)


purple_res <- goenrichment(colour="purple", allclusters = clustcolours, genenames = allGenes,
                           geneGOIDs = geneID2GO)


green_res <- goenrichment(colour="green", allclusters = clustcolours, genenames = allGenes,
                          geneGOIDs = geneID2GO)

lightcyan_res <- goenrichment(colour="lightcyan", allclusters = clustcolours, genenames = allGenes,
                              geneGOIDs = geneID2GO)

lightyellow_res <- goenrichment(colour="lightyellow", allclusters = clustcolours, genenames = allGenes,
                                geneGOIDs = geneID2GO)

darkslateblue_res <- goenrichment(colour="darkslateblue", allclusters = clustcolours, genenames = allGenes,
                                  geneGOIDs = geneID2GO)



coral1_res <- goenrichment(colour="coral1", allclusters = clustcolours, genenames = allGenes,
                           geneGOIDs = geneID2GO)

## example of how to summarize Revigo outfiles
revigo_summarize <- function(revigo, sig_go){
  revigo <- revigo %>% mutate(newdesc = ifelse(revigo$eliminated == 1, NA, revigo$description))
  revigo$newdesc <- na.locf(revigo$newdesc)
  rev_go <- cbind(revigo$term_ID, revigo$newdesc)
  colnames(rev_go) <- c("GO.ID", "Revigo")
  rev_merge <- merge(sig_go, rev_go,  by= "GO.ID") #%>% group_by(., newdesc) %>% 
  #    summarize(DEG=sum(Significant), Expected=sum(Expected))
  #  rev_merge <- rev_merge %>% melt(., id.vars="newdesc")
  return(rev_merge)
}

turq_rev <- read.csv("./data_files/turquoise_rev.txt")
turquoise_rev <- revigo_summarize(turq_rev, turquoise_res)

