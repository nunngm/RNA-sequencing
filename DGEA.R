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
# load PTI dataset
files <- file.path("jarad2020/counts", list.files("jarad2020/counts"))
samples <- c("hrc_h24_s1", "hrc_h24_s2", "hrc_h24_s3", "mock_h24_s1", "mock_h24_s2", "mock_h24_s3")
names(files) <- samples


##Sample information to build a model
infection = factor(c(rep("hrc", 3), rep("mock", 3)), levels = c("mock", "hrc")) #Treatment type, mg = mock solution, pst = P syringae pv tomato
hpi = factor(rep(24, 6), levels = c(24)) #Hours post infiltration


design_full <- data.frame(sample=names(files),
                          file=files,
                          infection =infection,
                          hpi=hpi,
                          batch = factor(c(1,2,2,1,2,2), levels = c(1,2))
)
design_full

##Very general model that just has all the sample information
model_full <- formula(~infection+batch)
rawData <- DESeqDataSetFromHTSeqCount(design_full,design=model_full)


## Grouping all the experimental variables into distinct treatment groups
rawData$group <- factor(paste0(rawData$infection,rawData$hpi, rawData$batch))
rawData$group = factor(rawData$group,levels=c("mock242","hrc242", "mock241", "hrc241"))

# Filter lowly expressed genes
keep <- rowMeans(counts(rawData)) >= 10 #Genes which on average have less than 10 reads
#keep = rowMeans(counts(rawData)[, 1:18]) >=10 &rowMeans(counts(rawData)[, 19:36]) >=10
rawData <- rawData[keep,]

# Build DESeq object based on distinct treatment groups
rawData@design = ~group
ptiData <- DESeq(rawData)

#Analyze the quality of the pti data set
hmcol = hcl_palettes(palette = "Berlin") #Setting the colour palatte
for_pca <- rlog(ptiData, blind=F)
rlogMat <- assay(for_pca) # just making a matrix of the counts that have been corrected for over-dispersion in a "blind" fashion

distsRL <- dist(t(rlogMat)) # Computes a distance matrix (Euclidian Distance)
mat <- as.matrix(distsRL)  # Make sure it is a matrix


infection = c(rep("hrcC-",3),rep("Mock",3))
infection <- as.factor(infection)
length(infection)
batch = factor(c(1,2,3,1,2,3), levels = c(1,2,3))

rownames(mat) <- colnames(mat) <-   with(colData(ptiData), paste( infection, paste0(hpi, "h"), batch, sep=":"))

hc <- hclust(distsRL
             ,method = "average"
)  # performs hierarchical clustering
par(mar=c(7,4,4,5)+0.1)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)# picking our colours

tiff(filename = "heatmap.tiff", height = 1080, width = 1400) #opens a tiff device
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=T, trace="none",
          cexRow = 2,
          cexCol = 2,
          col = rev(hmcol), margin=c(13,13)) #prints to tiff
dev.off()

## Plot the results (quickly but kinda nicely)
p =plotPCA(for_pca, ntop = 10000,
           intgroup=c("infection","hpi", "batch"))

p =p+geom_point(aes(size=1)) +guides(colour = guide_legend(override.aes = list(size = 8)))+theme(panel.grid.major = element_blank(), 
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
)

p
ggsave("myplot.pdf",plot = p) #output the PCA plot
pdf("silly.pdf")


##DGE analysis
sigGenes_pti= results(ptiData, contrast = c("group", "hrc24", "mock24"), tidy = T)
sigGenes_pti = sigGenes_pti[ is.na(sigGenes_pti$padj) == F,]
sigGenes_pti = sigGenes_pti[sigGenes_pti$log2FoldChange > 0.58  & sigGenes_pti$pvalue < 0.05,]
objectSymbol[sigGenes_pti$row]

sigGenes_m = sigGenes_m$row[sigGenes_m$log2FoldChange > 1]


#__________________________Mine2018 dataset________________
files <- file.path("mine2018/counts", list.files("mine2018/counts"))
samples <- c("dc3000_h12_1", "dc3000_h12_2", "dc3000_h12_3", "dc3000_h24_1", "dc3000_h24_2", "dc3000_h24_3", "mock_h12_1", "mock_h12_2", "mock_h12_3", "mock_h24_1", "mock_h24_2", "mock_h24_3", "rpm1_h12_1", "rpm1_h12_2", "rpm1_h12_3", "rpm1_h24_1", "rpm1_h24_2", "rpm1_h24_3", "rpt2_h12_1", "rpt2_h12_2", "rpt2_h12_3", "rpt2_h24_1", "rpt2_h24_2", "rpt2_h24_3")
names(files) <- samples


##Sample information to build a model
infection = factor(c(rep("dc3000", 6), rep("mock", 6), rep("rpm1", 6), rep("rpt2", 6)), levels = c("mock", "dc3000", "rpt2", "rpm1")) #Treatment type, mg = mock solution, pst = P syringae pv tomato
hpi = factor(rep(c(rep(12, 3), rep(24, 3)), 4), levels = c(12, 24)) #Hours post infiltration


design_full <- data.frame(sample=names(files),
                          file=files,
                          infection =infection,
                          hpi=hpi
)
design_full

##Very general model that just has all the sample information
model_full <- formula(~infection+hpi)
rawData <- DESeqDataSetFromHTSeqCount(design_full, design=model_full)


## Grouping all the experimental variables into distinct treatment groups
rawData$group <- factor(paste(rawData$infection, rawData$hpi, sep = "_"))
rawData$group = factor(rawData$group,levels=c("mock_12", "mock_24", "dc3000_12", "dc3000_24", "rpm1_12", "rpm1_24", "rpt2_12", "rpt2_24"))

# Filter lowly expressed genes
keep <- rowMeans(counts(rawData)) >= 10 #Genes which on average have less than 10 reads
#keep = rowMeans(counts(rawData)[, 1:18]) >=10 &rowMeans(counts(rawData)[, 19:36]) >=10
rawData <- rawData[keep,]

# Build DESeq object based on distinct treatment groups
rawData@design = ~group
etiData <- DESeq(rawData)

#Analyze the quality of the pti data set
hmcol = hcl_palettes(palette = "Berlin") #Setting the colour palatte
for_pca <- rlog(etiData, blind=T)
rlogMat <- assay(for_pca) # just making a matrix of the counts that have been corrected for over-dispersion in a "blind" fashion

distsRL <- dist(t(rlogMat)) # Computes a distance matrix (Euclidian Distance)
mat <- as.matrix(distsRL)  # Make sure it is a matrix

#Sample number labelling
batch = factor(rep(c(1,2,3), 8), levels = c(1,2,3))
rawData$batch = batch

rownames(mat) <- colnames(mat) <-   with(colData(etiData), paste( infection, paste0(hpi, "h"), sep=":"))

hc <- hclust(distsRL
             ,method = "average"
)  # performs hierarchical clustering
par(mar=c(7,4,4,5)+0.1)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)# picking our colours

tiff(filename = "heatmap.tiff", height = 1080, width = 1400) #opens a tiff device
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=T, trace="none",
          cexRow = 2,
          cexCol = 2,
          col = rev(hmcol), margin=c(13,13)) #prints to tiff
dev.off()

#----------------------------- Really good looking PCA for ETI data

for_pca <- rlog( etiData, blind = T )
rv <- rowVars(assay(for_pca))
# select the ntop genes by variance (across treatment groups)
ntop = 10000
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(for_pca)[select,]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
intgroup_df <- as.data.frame(colData(etiData)[, "group", drop = FALSE])


#Selecting the principle components
pc_x = 1
pc_y = 2
d <- data.frame(PC1 = pca$x[, pc_x], PC2 = pca$x[, pc_y], 
                group = intgroup_df, 
                inf = colData(for_pca)[,1],
                hpi = colData(for_pca)[,2],
                batch = colData(for_pca)[,4]#,
                #ageinf = as.integer(as.factor(paste0(d$age,d$)))
) #In pca$x[,]

#Setting the shapes for the various infection types
temp = c(rep(15,6), rep(3, 6), rep(16,6), rep(17,6))

#Drawing the PCA plot and demonstrating variance

#prints to tiff

ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "hpi")) + geom_point(size = 4, shape =temp, stroke = 1.5) +  xlab(paste0("PC ",pc_x," (", round(percentVar[pc_x] * 100), "% variance)")) + ylab(paste0("PC ",pc_y," (", round(percentVar[pc_y] * 100), "% variance)")) + 
        coord_fixed() +theme(panel.grid.major = element_blank(), 
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
        )



as.factor(paste0(d$inf, d$hpi))
transparency = c(rep(1,3), rep(1,3),rep(1,3), 
                 rep(1,3), rep(1,3),rep(1,3), 
                 rep(1,3), rep(1,3),rep(1,3), 
                 rep(1,3), rep(1,3),rep(1,3))
transparency = transparency[1:24]

tiff(filename = "PCA.tiff", height = 2000, width = 2000) #opens a tiff device
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "hpi")) + geom_point(size = 20,shape =temp, stroke = 6, alpha = transparency) +  xlab("") + ylab("") + coord_fixed() +theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                      panel.grid.minor = element_blank(),
                                                                                                                                                                                      panel.background = element_blank(), 
                                                                                                                                                                                      axis.line = element_line(colour = "black", size=4),
                                                                                                                                                                                      axis.title.x=element_text(size=15),
                                                                                                                                                                                      #axis.text.x=element_blank()),
                                                                                                                                                                                      axis.ticks=element_line(colour = "black", size =4),
                                                                                                                                                                                      axis.ticks.length = unit(20,"points") ,
                                                                                                                                                                                      axis.title.y = element_text(size=15),
                                                                                                                                                                                      legend.position = "right",
                                                                                                                                                                                      axis.text = element_text(size= 75),
)
dev.off()

ggsave("makegraph_noPoints.png")
pdf("myplots.pdf")

##Comparing ETI data to ARR data
sigGenes_dc3000 = results(etiData, contrast = c("group", "dc3000_12", "mock_12"),  alpha = 0.05, tidy = T)
sigGenes_dc3000 = sigGenes_dc3000[ is.na(sigGenes_dc3000$padj) == F,]
sigGenes_dc3000 = sigGenes_dc3000[abs(sigGenes_dc3000$log2FoldChange) > log2(1.5)  & sigGenes_dc3000$padj < 0.05,]

sigGenes_y = results(allData, contrast = c("group", "ypst12", "ymg12"), alpha = 0.05, tidy = T)
sigGenes_y = sigGenes_y[ is.na(sigGenes_y$padj) == F,]
sigGenes_y = sigGenes_y[sigGenes_y$log2FoldChange > 1  & sigGenes_y$padj < 0.05,]

sigGenes_m = results(allData, contrast = c("group", "mpst24", "mmg24"), alpha = 0.05, tidy = T)
sigGenes_m = sigGenes_m[ is.na(sigGenes_m$padj) == F,]
sigGenes_m = sigGenes_m[sigGenes_m$log2FoldChange > 1  & sigGenes_m$padj < 0.05,]

sigGenes_mmock = results(allData, contrast = c("group", "mmg12", "ymg12"), alpha = 0.05, tidy = T)
sigGenes_mmock = sigGenes_mmock[ is.na(sigGenes_mmock$padj) == F,]
sigGenes_mmock = sigGenes_mmock[sigGenes_mmock$log2FoldChange > 1  & sigGenes_mmock$padj < 0.05,]

sigGenes_eti = results(etiData, contrast = c("group", "rpm1_24", "mock_24"),  alpha = 0.05, tidy = T)
sigGenes_eti = sigGenes_eti[ is.na(sigGenes_eti$padj) == F,]
sigGenes_eti = sigGenes_eti[sigGenes_eti$log2FoldChange > 1  & sigGenes_eti$padj < 0.05,]

sigGenes_pti = results(ptiData, contrast = c("group", "flg_24hpt", "ctrl_24hpt"),  alpha = 0.05, tidy = T)
sigGenes_pti = sigGenes_pti[ is.na(sigGenes_pti$padj) == F,]
sigGenes_pti = sigGenes_pti[sigGenes_pti$log2FoldChange > 1  & sigGenes_pti$padj < 0.05,]

sigGenes_m = find.volcano("12", lfc.cut = 1)
sigGenes_m = sigGenes_m[sigGenes_m$de == "UP",]

sigGenes_y = find.volcano("12", lfc.cut = 1, SUGs = T)
sigGenes_y = sigGenes_y[sigGenes_y$de == "UP",]

colors = qualitative_hcl(4, palette = "dark",c = 90)
show_col(colors[1:4] )

venn.diagram(x = list(sigGenes_m$row, sigGenes_eti$row, sigGenes_pti$row),
             category.names = c("ARR", "ETI", "PTI"),
             filename = "12h ARR-ETI-PTI12hpt.tiff",
             output = T,
             imagetype = "tiff",
             euler.d = F,
             scaled = F,
             col =  colors[1:3],
             fill = "white",
             cat.col = "black",
             cat.cex = 2,
             cat.dist = c(0.12, 0.12,0.12),
             margin =0.15)

GOI = sigGenes_eti$row[ sigGenes_eti$row %in% sigGenes_pti$row]

GOI = rownames(sigGenes_m)[ rownames(sigGenes_m) %in% sigGenes_pti$row]

GOI = GOI[ GOI %in% sigGenes_pti$row]
GOI = GOI[ GOI %in% sigGenes_m$row]


#------ Mock vs Mock
rawData$group <- factor(paste0(rawData$age,rawData$infection))
rawData$group = factor(rawData$group,levels=c("ymg", "ypst", "mmg", "mpst"))

keep <- rowMeans(counts(rawData)) >= 10 #Genes which on average have less than 10 reads
#keep = rowMeans(counts(rawData)[, 1:18]) >=10 &rowMeans(counts(rawData)[, 19:36]) >=10
rawData <- rawData[keep,]

# Build DESeq object based on distinct treatment groups
rawData@design = ~group
arrData <- DESeq(rawData)

ageGenes = list(results(allData, contrast = c("group", "mmg0", "ymg0"), alpha = 0.05, tidy = T), results(allData, contrast = c("group", "mmg12", "ymg12"), alpha = 0.05, tidy = T), results(allData, contrast = c("group", "mmg24", "ymg24"), alpha = 0.05, tidy = T))

ageGenes = results(arrData, contrast = c("group", "mmg", "ymg"), alpha = 0.05, tidy = T)
ageGenes = ageGenes[ ageGenes$log2FoldChange > 1 & ageGenes$padj < 0.05 , ]


allComp = compare.group()
thresh = 0.05
lfc = 1

mpstcomp = list(results(allData,contrast = c("group", "mpst0", "ymg0"), alpha = 0.05, pAdjustMethod="BH", tidy = T), results(allData,contrast = c("group", "mpst12", "ymg12"), alpha = 0.05, pAdjustMethod="BH", tidy = T), results(allData,contrast = c("group", "mpst24", "ymg24"), alpha = 0.05, pAdjustMethod="BH", tidy = T))

#RLP protein fanily

mydata= read.table(file= "clipboard",sep= "\t",header =F)
mydata[, 2] = toupper(mydata[,2])
genes_of_interest = mydata[,2]
GOI = rownames(allComp$ypst0ymg0) %in% genes_of_interest

finalLFC = cbind(allComp$ypst0ymg0[GOI, 2], allComp$mmg0ymg0[GOI,2], mpstcomp[[1]][GOI, 3], allComp$ypst12ymg12[GOI, 2], allComp$mmg12ymg12[GOI, 2], mpstcomp[[2]][GOI, 3], allComp$ypst24ymg24[GOI, 2], allComp$mmg24ymg24[GOI, 2], mpstcomp[[3]][GOI, 3])

rownames(finalLFC) = objectSymbol[rownames(allComp$ypst0ymg0)[rownames(allComp$ypst0ymg0) %in% genes_of_interest]]
colnames(finalLFC) = c("Y.Pst", "M.Mock", "M.Pst","Y.Pst", "M.Mock", "M.Pst", "Y.Pst", "M.Mock", "M.Pst" )


#worked
hmcol1 = sequential_hcl(n = 10, h = 360, c =180, l = c(20,130), power = 0.4, rev = T)
show_col(hmcol1)
hmcol2 = sequential_hcl(n = 10, h = 230, c = 100, l = c(50,100), power = 0.6)
show_col(hmcol2)
hmcol = c(hmcol2[1:9], "#FFFFFF", hmcol1[2:10]) 
show_col(hmcol)


tiff("rplot2.tiff", width = 4000, height = 2000)
aheatmap(finalLFC, color = hmcol, Rowv = T, Colv = NA, distfun= "euclidean",treeheight = 200,hclustfun = "average", scale = "none", cellwidth = 320, breaks =0)
dev.off()

pdf()
heatmap.2(finalLFC, Rowv=T, Colv = F,
          trace="none", col = hmcol,
          dendrogram = "row", cexRow = 0.5)
dev.off()

