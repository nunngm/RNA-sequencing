# data = res$mpst24mmg24
# data$log2FoldChange = res$mmg24ymg24$log2FoldChange+res$mpst24mmg24$log2FoldChange-res$ypst24ymg24$log2FoldChange
# data$pvalue = (1-res$ypst24ymg24$padj)
res = compare.group(hpi= "12")

data = res[[4]]

data$padj =  -2*(log(1-res[[4]]$padj)+log(1-res[[2]]$padj)+log(res[[1]]$padj))
data$padj =  pchisq(data$padj,2*3) #Uses the inverse of fisher's method to basically discount samples which aren't: 1. Differentiall expressed in mature plants 2. differentially expressed in M.Pst compared to Y.Pst 3. Are not differentially expressed in Y.Pst compared to Y.Mock
data$log2FoldChange = res[[4]]$log2FoldChange

data$de = "NO"
data$de[data$log2FoldChange > 0.585 & data$padj < 0.05] <- "UP"
data$de[data$log2FoldChange < -0.585 & data$padj < 0.05] = "DOWN"
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
data$gene_symbol = getGeneName(rownames(data))
data$gene_symbol[is.na(data$gene_symbol)] = ""
data = data[is.na(data$log2FoldChange)==F & is.na(data$padj)==F,]

data$delabel = NA
n= 30
data$delabel[order(data$log2FoldChange * (data$padj <0.05) *(-log10(data$padj)), decreasing = T)][c(1:n,(nrow(data)-(n-1)):nrow(data))] = data$gene_symbol[order(data$log2FoldChange * (data$padj <0.05) *(-log10(data$padj)), decreasing = T)][c(1:n,(nrow(data)-(n-1)):nrow(data))]

data = res2

x = 10
nrow(tmp[[x]][is.na(tmp[[x]]$padj) == F & tmp[[x]]$padj<0.01 & tmp[[x]]$log2FoldChange>0, ])
nrow(tmp[[x]][is.na(tmp[[x]]$padj) == F & tmp[[x]]$padj<0.01 & tmp[[x]]$log2FoldChange<0, ])
# [data$de != "NO",][c(1:25, (nrow(data)-24):nrow(data)),]
# 
# temp = data[order(data$log2FoldChange, decreasing = T),]
# temp = temp[temp$padj < 0.05,]
# temp = temp[c(1:25, (nrow(temp)-24):nrow(temp)),]
# 
# data$delabel[rownames()] = data$gene_symbol
# data$delabel = 
# 
# 
# data$delabel[data$de != "NO"] = data$gene_symbol[data$de != "NO"]
data = find.volcano("24")
hour = c("0", "12", "24")
tmp = lapply(hour, find.volcano, de =1)

for (i in 1:length(tmp)){
  
}

write.xlsx( tmp[[1]], "volcano.xlsx", sheetName = "0 hpi",col.names = T, row.names = T, append = F)
write.xlsx( tmp[[2]], "volcano.xlsx", sheetName = "12 hpi",col.names = T, row.names = T, append = T)
write.xlsx( tmp[[3]], "volcano.xlsx", sheetName = "24 hpi",col.names = T, row.names = T, append = T)


tiff(filename = "Volcano_24h.tiff", height = 2000, width = 2000)
ggplot(data = data, aes(x = log2FoldChange, y = -log10(padj), color = factor(de, levels = c("UP", "NO", "DOWN")))) + geom_point(size = 5) +
  #geom_text_repel(col = "black") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black", size=2),
          axis.title.x=element_text(size=30),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =2), axis.ticks.length = unit(5,"points"),
          axis.title.y = element_text(size=30), legend.position = "right",
          axis.text = element_text(size=50)) +
    geom_vline(xintercept=c(-0.585, 0.585), col="black", linetype = 2, size = 2) + geom_hline(yintercept=-log10(0.05), col="black", linetype = 2, size = 2)+ 
          scale_color_manual(values=c("#DB282E", "#BEBEBE", "#2E70CC"))+ ylim(0,10)
dev.off()

library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=data, aes(x=log2FoldChange, y=-log10(padj), col=de, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
