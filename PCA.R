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

ggsave("myplot.pdf",plot = p) #output the PCA plot
pdf("silly.pdf")



## Plotting other PCA components
install.packages("ggrepel") #only need to run once
#This code was originally taken from the plotPCA function of DESeq2 but has been modifed so that any of the principal components can be viewed

library(ggrepel)

#Additional PCA components
for_pca <- rlog( allData, blind = T )
rv <- rowVars(assay(for_pca))
# select the ntop genes by variance (across treatment groups)
ntop = 10000
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(for_pca)[select,]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
intgroup_df <- as.data.frame(colData(allData)[, "group", drop = FALSE])
group <- if (length("group") > 1) {
  factor(apply(intgroup_df, 1, paste, collapse = " : "))
} else{
  colData(allData)[["group"]]
}

#Selecting the principle components
pc_x = 3
pc_y = 4
d <- data.frame(PC1 = pca$x[, pc_x], PC2 = pca$x[, pc_y], 
                group = intgroup_df, 
                age = colData(for_pca)[,1],
                inf = colData(for_pca)[,2],
                hpi = colData(for_pca)[,3],
                ageinf = as.integer(as.factor(paste0(d$age,d$)))) #In pca$x[,]
temp = as.integer(as.factor(paste0(d$inf,d$age)))
temp = lapply(temp, function(x){tmp = if (x>2){
    x=15+(x-3)
} else{
    x = x-1
  }
  return(tmp)
})
temp =unlist(temp)

#Drawing the PCA plot and demonstrating variance
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "hpi")) + geom_point(size = 4,shape =temp,stroke = 1.5) + xlab(paste0("PC ",pc_x," (", round(percentVar[pc_x] * 100), "% variance)")) + ylab(paste0("PC ",pc_y," (", round(percentVar[pc_y] * 100), "% variance)")) + coord_fixed() +theme(panel.grid.major = element_blank(), 
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

