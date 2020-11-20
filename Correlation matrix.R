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
# somethings depend on deseq2_drought
# uncomment if you want to source this file
#source("deseq2_drought.R)


allowWGCNAThreads()
options(stringsAsFactors = FALSE)

gsg <- goodSamplesGenes(t(allgenes), verbose = 3);
gsg$allOK
allgenes <- t(allgenes)
sampleTree = hclust(dist(allgenes), method = "average")

## Plot the sample tree: Open a graphic output window of size 12 by 9 inches
## The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
##pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#     cex.axis = 1.5, cex.main = 2)
#
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

## 8 looks like a good fit

# need to do it block wise because my computer doesnt have enough RAM for this
network = blockwiseModules(t(allgenes), power = 8,
                           TOMType = "none",
                           networkType = "signed hybrid",
                           randomSeed = 319, corType = "pearson", 
                           minModuleSize = 30,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "droughtTOM",
                           verbose = 3, maxBlockSize = 10000)

## if you save this, it is much quicker to run
#save.image("blockwiseNetwork.RData")
load("~/Documents/Manuscripts/drought/scripts/blockwiseNetwork.RData")


clustcolours <- labels2colors(network$colors)
MEs = moduleEigengenes(t(allgenes), clustcolours)$eigengenes




nGenes = nrow(allgenes)
nSamples = ncol(allgenes);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(t(allgenes), clustcolours)$eigengenes
MEs = orderMEs(MEs0)
ecotype_condition <- c(rep("YWW1", 4), 
                       rep("YD1", 3),
                       rep("YWW2", 4),
                       rep("YD2", 4),
                       rep("SWW1", 4), 
                       rep("SD1", 3),
                       rep("SWW2", 4),
                       rep("SD2", 5))

datTraits <- data.frame(ecocond = as.factor(ecotype_condition))
rownames(datTraits) <- colnames(allgenes)
traits <- model.matrix(~ ., data=datTraits,
                       contrasts.arg = lapply(datTraits, contrasts, contrasts=FALSE))
traits <- traits[,-1]
moduleTraitCor <- WGCNA::cor(MEs, traits, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,31);

pal <- colorRampPalette(c("#E7B800","#E7B800", "white", "#67A7C1", "#67A7C1"),
                        alpha=TRUE, interp="spline")(100)
droughtcol <- c("#C0C8C8","#8D9B9B",  "#6A90D1", "#5776AC",
                "#C0C8C8","#8D9B9B",  "#6A90D1", "#5776AC")

colnames(moduleTraitCor) <- c("SD1", "SD2", "SWW1", "SWW2", "YD1", "YD2", "YWW1", "YWW2")
sample_clust <- hclust(dist(t(moduleTraitCor), method="euclidean"))
module_clust <- hclust(as.dist(1-cor(t(moduleTraitCor), method="pearson")))

rownames(moduleTraitCor) <- rownames(moduleTraitCor) %>% substr(.,3,nchar(rownames(moduleTraitCor)))

head(moduleTraitCor)

colnames(moduleTraitCor) <- c("SD1", "SD2", "SWW1", "SWW2", "YD1", "YD2", "YWW1", "YWW2")

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


clust.of.interest <- c('lightcyan1', 'lightyellow', 'purple', 'pink', 
                       'blue', 'darkslateblue', 'turquoise', 'coral', 'red')

(heatmap <- ggplot(data = correlation_ggplot %>% filter(., Var1 %in% clust.of.interest), 
                   aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill =value), colour="darkgrey") +
    scale_fill_gradient2(midpoint = 0, low = "#8AA399", mid = "white", 
                         high = "#FF8977", na.value = "white", name = "Correlation", 
                         limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) + 
    labs(x="Ecotype and Condition", y="Clusters") +  theme_bw() + #theme(axis.text.y=element_blank(), plot.margin = margin(-0.75, 0, 0,0 , "cm")))
    theme(axis.text.y=element_text(angle=20,vjust=0.5), plot.margin = margin(-0.75, 0, 0,0 , "cm")))
(dendro_heatmap <- px + heatmap + plot_layout(ncol = 1, heights = c(1, 5))) # + plot_layout(ncol=2, widths=(1,5))

ggsave("./figs/select_clust_heatmap.pdf", dendro_heatmap,
       height=6, width = 7, units = "in")


(heatmap_full <- ggplot(data = correlation_ggplot, 
                        aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill =value), colour="darkgrey") +
    scale_fill_gradient2(midpoint = 0, low = "#00798C", #mid = "white", 
                         high = "#FF6F59", na.value = "white", name = "Correlation", 
                         limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) + 
    labs(x="Ecotype and Condition", y="Clusters") +  theme_bw() + #theme(axis.text.y=element_blank(), plot.margin = margin(-0.75, 0, 0,0 , "cm")))
    theme(axis.text.y=element_text(angle=20,vjust=0.5), plot.margin = margin(-0.75, 0, 0,0 , "cm")))
(dendro_heatmap_full <- px + heatmap_full + plot_layout(ncol = 1, heights = c(1, 5))) 


ggsave("./figs/all_clust_heatmap.pdf", dendro_heatmap_full,
       height=10, width = 8, units = "in")
geneClusters <- network$colors
genenames <- rownames(allgenes)
names(clustcolours) <- genenames

clustcolours[names(clustcolours) %>% toupper() %in% pattern3genes]
clustcolours[names(clustcolours) %>% toupper() %in% pattern4genes]


## What about ALL DEGS?? 
## let's put all DEGs into a list, then are able to see which cluster(s) entire groups of degs are in!
## this depends on output from deseq2_drought.R

SWW1SD1_alldegs <- rbind(SWW1.SD1.degs[['pos']], SWW1.SD1.degs[['neg']])
SD1SWW2_alldegs <- rbind(SD1.SWW2.degs[['pos']], SD1.SWW2.degs[['neg']])
SWW2SD2_alldegs <- rbind(SWW2.SD2.degs[['pos']], SWW2.SD2.degs[['neg']])
YWW1YD1_alldegs <- rbind(YWW1.YD1.degs[['pos']], YWW1.YD1.degs[['neg']])
YD1YWW2_alldegs <- rbind(YD1.YWW2.degs[['pos']], YD1.YWW2.degs[['neg']])
YWW2YD2_alldegs <- rbind(YWW2.YD2.degs[['pos']], YWW2.YD2.degs[['neg']])

alldegs_list <- list(sww1sd1 = SWW1SD1_alldegs, sd1sww2=SD1SWW2_alldegs, sww2sd2=SWW2SD2_alldegs,
                     yww1yd1 = YWW1YD1_alldegs, yd1yww2=YD1YWW2_alldegs, yww2yd2=YWW2YD2_alldegs)



test<-SWW1SD1_alldegs$Row.names
clustcolours[names(clustcolours) %in% test]


clust_count <- clustcolours[names(clustcolours) %in% alldegs] %>% 
  table() # %>% as.matrix() %>% sort() %>% rev()

clust_degs <- clustcolours[names(clustcolours) %in% alldegs] %>% table() %>% as.matrix()
clustnum <- clustcolours %>% table() %>% as.matrix
perc_deg <- merge(clust_degs, clustnum, all.y=T, by=0)
perc_deg <- data.frame(perc_deg, perc=(perc_deg$V1.x/perc_deg$V1.y)*100)
perc_deg <- merge(perc_deg, moduleTraitCor, by.x="Row.names",by.y=0)%>% arrange(.,desc(perc))

#write.csv(perc_deg, "deg_clusts.csv")
## how many clusters do we have to do to get 90% of the degs?
clust_count[1:18] %>% sum() ## too many..what if you were just doing the ecotypes ind?

y.clust_count <- clustcolours[names(clustcolours) %in% yukdeg] %>% 
  table() %>% as.matrix() %>% sort() %>% rev() ## 3559

y.clust_count[1:6] %>% sum()

s.clust_count <- clustcolours[names(clustcolours) %in% shandeg] %>% 
  table() %>% as.matrix() %>% sort() %>% rev() ##5580
s.clust_count[1:5] %>% sum

## does this change with lncRNAs?

SWW1SD1lnc <- rbind(SWW1.SD1.degs[['pos_lnc']], SWW1.SD1.degs[['neg_lnc']])
SD1SWW2lnc <- rbind(SD1.SWW2.degs[['pos_lnc']], SD1.SWW2.degs[['neg_lnc']])
SWW2SD2lnc <- rbind(SWW2.SD2.degs[['pos_lnc']], SWW2.SD2.degs[['neg_lnc']])
YWW1YD1lnc <- rbind(YWW1.YD1.degs[['pos_lnc']], YWW1.YD1.degs[['neg_lnc']])
YD1YWW2lnc <- rbind(YD1.YWW2.degs[['pos_lnc']], YD1.YWW2.degs[['neg_lnc']])
YWW2YD2lnc <- rbind(YWW2.YD2.degs[['pos_lnc']], YWW2.YD2.degs[['neg_lnc']])

lncdegs <- list(SWW1SD1=SWW1SD1lnc, SD1SWW2=SD1SWW2lnc, SWW2SD2=SWW2SD2lnc,
                YWW1YD1=YWW1YD1lnc, YD1YWW2=YD1YWW2lnc, YWW2YD2=YWW2YD2lnc)

lncdegs_genes <- rbind(SWW1SD1lnc, SD1SWW2lnc,SWW2SD2lnc,
                       YWW1YD1lnc,YD1YWW2lnc,YWW2YD2lnc)
lncdegs_genes <- lncdegs_genes$Row.names %>% unique()
lnc_clusters <- lapply(lncdegs, function(x) clustcolours[names(clustcolours) %in% x$Row.names])

lnc_clusters[["SD1SWW2"]] %>% table() %>% sort() 
lnc_clusters[["SD1SWW2"]] %>% table() %>% sort()
lnc_clusters[["SWW2SD2"]] %>% table() %>% sort()             
lnc_clusters[["YD1YWW2"]] %>% table() %>% sort() 
lnc_clusters[["YD1YWW2"]] %>% table() %>% sort()
lnc_clusters[["YWW2YD2"]] %>% table() %>% sort()             



lnc_count <- clustcolours[names(clustcolours) %in% lncdegs_genes] %>% 
  table() # %>% as.matrix() %>% sort() %>% rev()




lnc_degs <- clustcolours[names(clustcolours) %in% lncdegs_genes] %>% table() %>% as.matrix()
clustnum <- clustcolours %>% table() %>% as.matrix
perc_lnc <- merge(lnc_degs, clustnum, all.y=T, by=0)
perc_lnc <- data.frame(perc_lnc, perc=(perc_lnc$V1.x/perc_lnc$V1.y)*100)
perc_lnc <- merge(perc_lnc, moduleTraitCor, by.x="Row.names",by.y=0)%>% arrange(.,desc(perc))

# how many genes in each cluster are lncRNAs?
total_lnc_count <- clustcolours[names(clustcolours) %in% unique_lncrna] %>% table() %>% as.matrix()
clustnum <- clustcolours %>% table() %>% as.matrix()
total_perc_lnc <- merge(total_lnc_count, clustnum, by=0)
total_perc_lnc <- data.frame(total_perc_lnc, perc=(total_perc_lnc$V1.x/total_perc_lnc$V1.y)*100)
# how many genes in each cluster are not in reference annotation?
novel_genes <-names(clustcolours)[!startsWith(names(clustcolours), "Thhalv")]
novel_count <- clustcolours[names(clustcolours) %in% novel_genes] %>% table() %>% as.matrix()
novel_perc <- merge(novel_count, clustnum, by=0)
novel_perc <- data.frame(novel_perc, perc=(novel_perc$V1.x/novel_perc$V1.y)*100)

clustcolours[grep("Thhalv10012372", names(clustcolours))]
clustcolours[grep("Thhalv10023284", names(clustcolours))]

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

## DRIR analysis
## are any of our drought related DEGs homologous to DRIR (the drought related lncRNA in ath?)

drir_blast <- read.csv("./data_files/yuk_drir.out",
                       sep="\t", header=F)
gene_trans <- read.csv("./data_files/gene_transcript.csv",
                       header=F, sep=" ")

drir_trans <- merge(drir_blast, gene_trans, by.x=2, by.y=2, all.x=T)

drir_trans <- drir_trans[!duplicated(drir_trans$V1.y),]

drir_trans$V1.y %in% YWW1.YD1.degs[["pos"]]$Row.names
drir_trans$V1.y %in% YWW2.YD2.degs[["pos"]]$Row.names
drir_trans$V1.y %in% SWW1.SD1.degs[["pos"]]$Row.names
SWW2.SD2.degs[["pos"]][SWW2.SD2.degs[["pos"]]$Row.names %in% drir_trans$V1.y,] 

turquoise_genes[turquoise_genes %in% drir_trans$V1.y]

a <- turquoise_genes[turquoise_genes %in% YWW1.YD1.degs[["pos_lnc"]]$Row.names]
b <- turquoise_genes[turquoise_genes %in% YWW2.YD2.degs[["pos_lnc"]]$Row.names]
turquoise_genes[turquoise_genes %in% SWW1.SD1.degs[["pos_lnc"]]$Row.names]
c <- turquoise_genes[turquoise_genes %in% SWW2.SD2.degs[["pos_lnc"]]$Row.names]

c(a,b,c) %>% unique()
yuk_droughtlnc <- c(a,b)

intersect(yuk_droughtlnc, c)
setdiff(yuk_droughtlnc, c)
setdiff(c, yuk_droughtlnc)