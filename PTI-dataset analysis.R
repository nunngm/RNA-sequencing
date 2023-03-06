

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
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Proj-PTI RNA-seq\\Exp-R workshop")
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-PTI RNA-seq\\Exp-R workshop")
#Desktop directory
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Proj-PTI RNA-seq\\Exp-R workshop")


gene_associations <- read.delim("gene_association_final.txt", comment.char = "!", header = FALSE, as.is = TRUE) 
colnames(gene_associations) <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
                                 "DB:Reference", "Evidence", "With_From", "Aspect", "DB_Object_Name",
                                 "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_by")

# # Uncomment if this is the first time running this code this helps produce a gene_association file to work with in the future
# #I didn't want to look at locus ID's from TAIR so I took the first item of DB_Object_Synonym which is the gene ID and put it in the DB_Object_ID column as that is more useful
# gene_associations$DB_Object_ID = res
# gen_association_save <- sapply(gene_associations, function(x){paste(x, collapse = ", ")}) 
# gen_association_save <- data.frame("Gene Association" = gene_associations)
# write.table(gen_association_save, file = "gene_association_final.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = F)

# Trimming the dataframe so that it's only what we're interested in
gene_associations <- gene_associations[,c(2,3,5,7,9,14)]


##Make Gene-Go object from scratch
# # Go through every unique gene and pull out any GO_ID associated with the given gene then give
# gene_GO <- lapply(unique(gene_associations$DB_Object_ID), function(x){tmp <- gene_associations %>% filter(DB_Object_ID == x)
# return(tmp$GO_ID)})
# names(gene_GO) <- unique(gene_associations$DB_Object_ID) 
# #Save it
# gene_GO_save <- sapply(gene_GO, function(x){paste(x, collapse = ", ")}) 
# gene_GO_save <- as.data.frame( gene_GO_save) # Making a dataframe
# write.table(gene_GO_save, file = "TAIR_to_GO.delim", sep = "\t", quote = FALSE, col.names = FALSE) 


#Load a pre-made Gene-to-go from file
gene_GO <- readMappings("TAIR_to_GO.delim")

#Collect filenames and label samples
files <- file.path("counts", list.files("counts"))
samples <- as.character(as.data.frame(strsplit(list.files("counts"), split = ".", fixed = T))[1,]) #split the string name from the rest of the files
names(files) <- samples

##Sample information to build a model
treatment = factor(c(rep("flg", 24), rep("mock", 24), rep("unt", 3)), levels = c("unt", "mock", "flg")) #Treatment type, unt = untreated, mock = mock treatment (sterile H2O), flg = treated with flg22

#Due to untreated being nested within itself changing mock and unt into one group called control
treatment = factor(c(rep("flg", 24), rep("ctrl", 24), rep("ctrl", 3)), levels = c("ctrl", "flg"))

hpi = factor(c(rep(c(rep("12hpi", 3), rep("12hpt",3), rep("18hpi", 3), rep("18hpt",3), rep("24hpi", 3), rep("24hpt", 3), rep("6hpi", 3), rep("6hpt", 3)) , 2) ,rep("0hpt", 3)), levels = c("0hpt", "6hpt", "12hpt", "18hpt", "24hpt", "6hpi", "12hpi", "18hpi", "24hpi")) # timepoint entry: hpt = hours post treatment, hpi = hours post infection. Plants were treated with treatment denoted above then after 24 hours were infected with Pst DC3000 pVSP61

design_full <- data.frame(sample=names(files),
                          file=files,
                          treatment = treatment,
                          hpi=hpi
)
design_full

##Very general model that just has all the sample information
model_full <- formula(~treatment+hpi)
rawData <- DESeqDataSetFromHTSeqCount(design_full,design=model_full)


## Grouping all the experimental variables into distinct treatment groups
rawData$group <- factor(paste(rawData$treatment,rawData$hpi, sep = "_"))
rawData$group = factor(rawData$group,levels=c("ctrl_0hpt","ctrl_6hpt", "flg_6hpt", "ctrl_12hpt", "flg_12hpt", "ctrl_18hpt", "flg_18hpt", "ctrl_24hpt", "flg_24hpt", "ctrl_6hpi", "flg_6hpi", "ctrl_12hpi", "flg_12hpi", "ctrl_18hpi", "flg_18hpi", "ctrl_24hpi", "flg_24hpi"))

# Filter lowly expressed genes
keep <- rowMeans(counts(rawData)) >= 10 #Genes which on average have less than 10 reads
#keep = rowMeans(counts(rawData)[, 1:18]) >=10 &rowMeans(counts(rawData)[, 19:36]) >=10
rawData <- rawData[keep,]

# Build DESeq object based on distinct treatment groups
rawData@design = ~group
ptiData <- DESeq(rawData)

##Rudimentary analysis of the dataset
hmcol = hcl_palettes(palette = "Berlin") #Setting the colour palatte
for_pca <- rlog(ptiData, blind=T)
rlogMat <- assay(for_pca) # just making a matrix of the counts that have been corrected for over-dispersion in a "blind" fashion

distsRL <- dist(t(rlogMat)) # Computes a distance matrix (Euclidian Distance)
mat <- as.matrix(distsRL)  # Make sure it is a matrix


rownames(mat) <- colnames(mat) <-   with(colData(etiData), paste( infection, paste0(hpi, "h"), sep=":"))

hc <- hclust(distsRL
             ,method = "average"
)  # performs hierarchical clustering
par(mar=c(7,4,4,5)+0.1)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)# picking our colours

tiff(filename = "heatmap.tiff", height = 2000, width = 2000) #opens a tiff device
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=T, trace="none",
          cexRow = 2,
          cexCol = 2,
          col = rev(hmcol), margin=c(13,13)) #prints to tiff
dev.off()

rv <- rowVars(assay(for_pca))
# select the ntop genes by variance (across treatment groups)
ntop = 10000
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(for_pca)[select,]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
intgroup_df <- as.data.frame(colData(ptiData)[, "group", drop = FALSE])
group <- if (length("group") > 1) {
  factor(apply(intgroup_df, 1, paste, collapse = " : "))
} else{
  colData(ptiData)[["group"]]
}

#Selecting the principle components
pc_x = 1
pc_y = 2
d <- data.frame(PC1 = pca$x[, pc_x], PC2 = pca$x[, pc_y], 
                group = intgroup_df, 
                treatment = colData(for_pca)[,1],
                hpi = colData(for_pca)[,2]#,
                #ageinf = as.integer(as.factor(paste0(d$age,d$)))
) #In pca$x[,]
temp = as.integer(as.factor(paste0(d$treatment,d$hpi)))
temp = c(rep(c(rep(16, 3), rep(1,3), rep(17, 3), rep(2,3), rep(18, 3), rep(5, 3), rep(15, 3), rep(0, 3)) , 2) ,rep(3, 3))


#Drawing the PCA plot and demonstrating variance

#prints to tiff


as.factor(paste0(d$inf,d$age, d$hpi))
transparency = c(rep(c(rep(1, 3), rep(1,3), rep(1, 3), rep(1,3), rep(1, 3), rep(1, 3), rep(1, 3), rep(1, 3)) , 2) ,rep(1, 3))


tiff(filename = "PCA.tiff", height = 2000, width = 2000) #opens a tiff device
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "treatment")) + geom_point(size = 15,shape =temp, stroke = 6, alpha = transparency) +  xlab("") + ylab("") + coord_fixed() +theme(panel.grid.major = element_blank(), 
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

view.pti.gene = function()
