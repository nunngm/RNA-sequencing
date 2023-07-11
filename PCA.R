#PCA component analysis

for_pca <- rlog(allData, blind=T) #Blind = F as we are trying to see if something is wrong with the experiment design so we want rlog to know all the variables to see if there is a crazy outlier
dim(for_pca) ##Making sure all genes and samples are preserved


##very clean looking PCA plot
par(cex=0.2)
p =plotPCA(for_pca, ntop = 10000,
           intgroup=c("age","infection","hpi"))

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



## Plotting other PCA components
install.packages("ggrepel") #only need to run once
#This code was originally taken from the plotPCA function of DESeq2 but has been modifed so that any of the principal components can be viewed

library(ggrepel)

#Additional PCA components


for_pca <- rlog( allData, blind = T )
sel_pca = for_pca[,sub_24]
rv <- rowVars(assay(sel_pca))



# select the ntop genes by variance (across treatment groups)
ntop = 1000
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

PCAplot = function(rlogDat, #regularized log full data set
                   genes = 1000, # can either be a single number, which selects the top x varying genes, a list of characters, which selects rownames which match items from the user entered list, or can be a list of numbers which will just select the rows as numbered
                   samples = c(1:36), #defaults to all samples
                   pc_x = 1, #the PC on the x-axis
                   pc_y = 2, #the PC on the y-axis,
                   width = 5, height = 4,
                   colours = c("red", "blue"),
                   shapes = c(1,16,0,15),
                   print2file = F,
                   colourVar = "age",
                   shapeVar = "ageinf",
                   legend.position = "none"
){
  rlogDat = rlogDat[,samples]
  
  if (length(genes) == 1){
    rv <- rowVars(assay(rlogDat))
    genes <- order(rv, decreasing=TRUE)[seq_len(min(genes, length(rv)))] #selects the 'genes' number of top varying genes, if genes>the number of rows then just defaults to all genes
  }else if(is.character(genes)) {
    genes = rownames(assay(rlogDat)) %in% genes
  }
  pca <- prcomp(t(assay(rlogDat)[genes,]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  intgroup_df <- as.data.frame(colData(rlogDat)[, "group", drop = FALSE])
  group <- if (length("group") > 1) {
    factor(apply(intgroup_df, 1, paste, collapse = " : "))
  } else{
    colData(rlogDat)[["group"]]
  }
  d <- data.frame(PC1 = pca$x[, pc_x], PC2 = pca$x[, pc_y], 
                  group = intgroup_df, 
                  age = colData(rlogDat)[,1],
                  inf = colData(rlogDat)[,2],#,
                  hpi = colData(rlogDat)[,3],
                  ageinf = factor(paste0(colData(rlogDat)[,1],colData(rlogDat)[,2]), levels = c("ymg", "ypst", "mmg", "mpst"))
  )
  
  p = ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = colourVar), ) + 
    geom_point(size = 4,shape = shapes[as.integer(unlist(d[shapeVar]))], stroke = 1.5) +  
    xlab(paste0("PC ",pc_x," (", round(percentVar[pc_x] * 100), "% variance)")) + 
    ylab(paste0("PC ",pc_y," (", round(percentVar[pc_y] * 100), "% variance)")) + 
    #coord_fixed() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black", size=1),
          axis.title.x=element_text(size=15),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length = unit(5,"points") ,
          axis.title.y = element_text(size=15),
          legend.position = legend.position,
          axis.text = element_text(size=15),
    ) + scale_color_manual(values = colours)
  if (print2file == T){
    fileName = readline(prompt = "Enter the name of the file: ")
    ggsave(file = paste0(fileName,".svg"), plot = p, width = width, height = height)
  }else{p}
}

genes = 

# perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(sel_pca)[select,]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  intgroup_df <- as.data.frame(colData(sel_pca)[, "group", drop = FALSE])
  group <- if (length("group") > 1) {
    factor(apply(intgroup_df, 1, paste, collapse = " : "))
  } else{
    colData(sel_pca)[["group"]]
  }
  
  

#Selecting the principle components
pc_x = 1
pc_y = 2
d <- data.frame(PC1 = pca$x[, pc_x], PC2 = pca$x[, pc_y], 
                group = intgroup_df, 
                age = colData(sel_pca)[,1],
                inf = colData(sel_pca)[,2]#,
                #hpi = colData(sel_pca)[,3]#,
                #ageinf = as.integer(as.factor(paste0(d$age,d$)))
     ) #In pca$x[,]
temp = as.integer(as.factor(paste0(d$inf,d$hpi)))
temp = lapply(temp, function(x){tmp = if (x>2){
    x=15+(x-3)
} else{
    x = x-1
  }
  return(tmp)
})
temp =unlist(temp)
temp = c(1,1,1,16,16,16,0,0,0,15,15,15)


#Drawing the PCA plot and demonstrating variance

 #prints to tiff

p = ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "age"), ) + geom_point(size = 4,shape =temp, stroke = 1.5) +  xlab(paste0("PC ",pc_x," (", round(percentVar[pc_x] * 100), "% variance)")) + ylab(paste0("PC ",pc_y," (", round(percentVar[pc_y] * 100), "% variance)")) + 
    #coord_fixed() +
  theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=1),
    axis.title.x=element_text(size=15),
    #axis.text.x=element_blank()),
    axis.ticks=element_line(colour = "black", size =1),
    axis.ticks.length = unit(5,"points") ,
    axis.title.y = element_text(size=15),
    legend.position = "none",
    axis.text = element_text(size=15),
) + scale_color_manual(values = c("#0993AE", "#54B031" ))
p
ggsave(file = "24hpi.svg", plot = p, width = 5, height = 4)


as.factor(paste0(d$inf,d$age, d$hpi))
transparency = c(rep(1,3), rep(1,3),rep(1,3), 
                 rep(1,3), rep(1,3),rep(1,3), 
                 rep(1,3), rep(1,3),rep(1,3), 
                 rep(1,3), rep(1,3),rep(1,3))
transparency = transparency[1:24]

tiff(filename = "PCA_microarray.tiff", height = 2000, width = 2000) #opens a tiff device
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
