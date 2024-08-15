# Graphing a bunch of grene transcrtipts
## Note: assumes you have run ARR-RNA setup script before beginning

## Required libraries
library(tidyr)
library(tidyverse)

# load in genes of interest
mydata= read.table(file= "clipboard",sep= "\t",header =T) # three columns labelled, label = gene name as you want it displaued, accession = the accession number, annotation = for later annotation (depreciated)

fullCount = counts(allData, normalized = T)
fullCount = fullCount[rownames(fullCount) %in% mydata$accession,]

cull = colnames(fullCount)
cull = as.data.frame(t(as.data.frame(strsplit(cull, split = "_"))))
colnames(cull) = c("age", "treatment", "hpi", "biorep")

fullCount = fullCount[rownames(fullCount) %in% mydata$accession, cull$hpi == "12h"]
fullCount = t(fullCount)
sampleGroups = as.data.frame(t(as.data.frame(strsplit(rownames(fullCount), split = "_"))))
colnames(sampleGroups) = c("age", "treatment", "hpi", "biorep")
sampleGroups$age = factor(sampleGroups$age, levels = c("y", "m"))
sampleGroups$treatment = factor(sampleGroups$treatment, levels = c("mg", "pst"))
finalCount = cbind(sampleGroups,fullCount)

write.table(finalCount, "clipboard", sep = "\t", row.names = T)


## Graphing
df.avg = read.table("clipboard", sep = "\t",header=T)
df.avg$ageTreat = factor(df.avg$ageTreat, levels = c("y_mg", "y_pst", "m_mg", "m_pst"))
df.avg$age = factor(df.avg$age, levels = c("y", "m"))
df.avg$treatment = factor(df.avg$treatment, levels = c("mg", "pst"))
ref = df.avg[,colnames(df.avg)=="AT1G76850"]
for (i in 6:ncol(df.avg)){
  df.avg[,i] = df.avg[,i]/ref
}
df.avg = df.avg[,colnames(df.avg)!="AT1G76850"]
df.long = gather(df.avg,gene, expression,-biorep, -treatment, -age, -hpi, -ageTreat)
df.long$expression[df.long$expression==0] = NA
df.long$gene = factor(df.long$gene, levels = mydata$accession)
# for (i in 1:nrow(df.long)){
#   df.long$expression[i] = df.long$expression[i]*mydata[mydata$accession==df.long[i,]$gene,]$transcriptLength/1000
# }

sampCol = c("#FFFFFF", "#FF3853", "#00BBFF", "#5D2169")

p = ggplot(df.long, aes(x=gene, y=expression, fill = ageTreat )) + stat_summary(fun = mean, position = position_dodge(width = 0.75), geom = "bar", colour = "#000000", size = 0.75, width = 0.65) +
  stat_summary(fun = mean,
               fun.min = function(x) {mean(x) - sd(x)},
               fun.max = function(x) {mean(x) + sd(x)},
               geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000", position = position_dodge(width = 0.75)) +
  geom_jitter( color= "#000000",
               size=2, alpha=0.5, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75)) +
  scale_x_discrete(labels = mydata$label) + scale_y_continuous(limits = c(0,23), expand = c(0,0)) +
  scale_fill_manual(values = sampCol) +
  xlab("") + ylab("Normalized transcript abundance (relative to SEC5A)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size=1),
        axis.title.x=element_text(size=15),
        axis.text.x=element_text(face = "italic"),
        axis.ticks.y=element_line(colour = "black", size =1),
        axis.ticks.length.y = unit(5, "points"),
        axis.ticks.length.x = unit(5, "points"),
        axis.ticks.x=element_line(colour = "black", size = 1),
        #ggh4x.axis.ticks.length.minor = rel(5),
        axis.ticks.length = unit(5,"points") ,
        axis.title.y = element_text(size=15),
        axis.text = element_text(color = "black", size=15),
        legend.position = "none"
  )
p

setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR-NPR&PEN3\\graphs")

ggsave(file = "PDRexpression-RNAseq.svg", plot = p, width = 14, height = 7)
