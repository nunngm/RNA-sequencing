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
mydata = gather(mydata, genotype, measurement, Y.Col:M.b1b1)
mydata$measurement = log10(mydata$measurement)
mydata = mydata[!is.na(mydata$measurement),]

mydata$genotype = factor(mydata$genotype, levels = c("Y.Col", "M.Col", "Y.ald1", "M.ald1",  "Y.sard4", "M.sard4", "Y.fmo1", "M.fmo1"))

mydata$genotype = factor(mydata$genotype, levels = c("Y.Col", "M.Col",  "Y.bak1_5", "M.bak1_5", "Y.bkk1_1", "M.bkk1_1", "Y.b1b1","M.b1b1"))

mydata$genotype = factor(mydata$genotype, levels = c("Y.Col", "M.Col",  "Y.pen3.4", "M.pen3.4", "Y.pdr12.3", "M.pdr12.3", "Y.p3p12","M.p3p12"))
# Plot
mydata %>%
  ggplot( aes(x=genotype, y=measurement )) +
    geom_boxplot(fill = rep(c("#595959","#FFFFFF"), 4)) +
    geom_jitter( width = 0.25, color=as.integer(as.factor(mydata$experiment)),
                  size=2, alpha=0.5)  + coord_cartesian(ylim = c(4, 7.2))+scale_y_continuous(breaks=c(4, 5, 6, 7), 
                   labels = expression(10^4, 10^5, 10^6, 10^7)) +   theme(
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
svg(filename = "CSRARR.svg", width = 5, height = 4)
dev.off()

make_labels = function(breaks){
     labels = bquote(paste0('10^',as.character(breaks[1]))
     return(labels)
}

mydata %>%
  ggplot( aes(x=genotype, y=measurement )) +
     stat_summary(fun = mean, geom = "bar", fill = rep(c("#595959","#FFFFFF"), 4), colour = "#000000", size = 0.75) +
     stat_summary(fun = mean,
          fun.min = function(x) {mean(x) - sd(x)}, 
          fun.max = function(x) {mean(x) + sd(x)}, 
          geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000") +
    geom_jitter( width = 0.25, color=as.integer(as.factor(mydata$experiment)),
                  size=2, alpha=0.5) + coord_cartesian(ylim = c(4, 7.5)) +    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) + 
     scale_y_continuous(breaks=c(4, 5, 6, 7), 
                   labels = expression(10^4, 10^5, 10^6, 10^7),
                   expand = c(0, 0))+
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

analysis = mydata[,2:3]
colnames(analysis) = c("type", "cfu")
