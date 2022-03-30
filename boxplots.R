## Making Boxplots

# Libraries
library(tidyverse)
library(tidyr)
library(ggplot2)

#laptop directory
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")

#Desktop directory
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")

mydata= read.table(file= "clipboard",sep= "\t",header =T)
mydata = gather(mydata, genotype, measurement, Y.Col:M.fmo1)
mydata$measurement = log10(mydata$measurement)
mydata = mydata[!is.na(mydata$measurement),]

mydata$genotype = factor(mydata$genotype, levels = c("Y.Col", "M.Col", "Y.ald1", "M.ald1",  "Y.sard4", "M.sard4", "Y.fmo1", "M.fmo1"))
# Plot
mydata %>%
  ggplot( aes(x=genotype, y=measurement )) +
    geom_boxplot(fill = rep(c("#595959","#FFFFFF"), 4)) +
    geom_jitter( width = 0.25,color=as.character(mydata$experiment), size=2, alpha=0.7) + ylim(4,7.5)+
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    xlab("Age and genotype") + ylab(bquote('Bacterial level (CFU leaf disc'^-1*')')) + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=1),
    axis.title.x=element_text(size=15),
    #axis.text.x=element_blank()),
    axis.ticks=element_line(colour = "black", size =1),
    axis.ticks.length = unit(5,"points") ,
    axis.title.y = element_text(size=15),
    axis.text = element_text(color = "black", size=15)
)

svg(filname = "PIPARR.svg", width = 7, height = 7)
dev.off()