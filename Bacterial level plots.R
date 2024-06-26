## Making Boxplots

# Libraries
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(grDevices)
library(agricolae)
library(dplyr)
detach(package:plyr)

#laptop directory
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR-NPR&PEN3\\exp-PEN3 ARR assays")
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR-NPR&PEN3\\exp-NPR ARR assay")
#Desktop directory
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")


YMBQGraph = function(data, #your data in long format
                     ageCol = c("#595959","#FFFFFF"), # colour of bars, defaults to dark grey and white
                     ylim = c(4,8), #y axis limits (10^x), by default 10^4, 10^8
                     expCol = NA, #if you have multiple experiments you can colour the different points different colours, here is where you assign the colours
                     graph = F, #should you output the graph to a file, by default it outputs to R, the output file is .svg which is a vectored image (meanaing it can be scaled up or down infinitely)
                     width = 5, height = 4, #size of the output plot
                     exptID = "temp", #for the output file name, put the id of the experiment
                     barLabs,
                     box = F #finally, graph defaults to a bar graph (mean+error bars), but can display a boxplot by setting this to True
                     ){
  data = data %>% mutate(sampGroup = paste(age, genotype, sep = "_"), cfu = log10(cfu), .keep = "all")
  if(is.na(expCol)){expCol = as.integer(as.factor(data$experiment))}
  if(length(barLabs)<2){barLabs = as.character(levels(data$genotype))}
  faces = c("plain", rep("italic", times = length(levels(data$genotype))-1))
  print(data)
  anovaModel = aov(cfu ~ sampGroup, data = data)
  data %>% group_by(genotype,age)%>% summarise(cfu = mean(cfu)) %>% summarise(foldDiff = 10^(cfu[age == 'Y'] - cfu[age == 'M']))
  print(HSD.test(anovaModel, alpha=0.05, "sampGroup", console=F)$groups)
  
  p = ggplot(data, aes(x=genotype, y=cfu, fill = age ))
  if(box == T ){
    p = p +geom_boxplot(position = position_dodge(width = 0.9), width = 0.8, size = 0.75)
    }else{
    p = p + stat_summary(fun = mean, position = position_dodge(width = 0.75), geom = "bar", colour = "#000000", size = 0.75, width = 0.65) +
      stat_summary(fun = mean,
                   fun.min = function(x) {mean(x) - sd(x)},
                   fun.max = function(x) {mean(x) + sd(x)},
                   geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000", position = position_dodge(width = 0.75))
    }
  p = p + geom_jitter( color= expCol,
                 size=2, alpha=0.5, position = position_jitterdodge(dodge.width = 0.75)) + coord_cartesian(ylim = ylim) +    theme(
                   legend.position="none",
                   plot.title = element_text(size=11)
                 ) + 
    scale_y_continuous(breaks=c(4, 5, 6, 7, 8), 
                       labels = expression(10^4, 10^5, 10^6, 10^7, 10^8),
                       expand = c(0, 0.05))   + scale_x_discrete(labels = barLabs) +
    scale_fill_manual(values = ageCol) +
    xlab("Genotype") + ylab(bquote('Bacterial level (cfu leaf disc'^-1*')')) + 
    theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black", size=1),
         axis.title.x=element_text(size=15),
         axis.text.x=element_text(face = faces),
         axis.ticks.y=element_line(colour = "black", size =1),
         axis.ticks.length.y = unit(5, "points"),
         axis.ticks.length.x = unit(5, "points"),
         axis.ticks.x=element_line(colour = "black", size = 1),
         #ggh4x.axis.ticks.length.minor = rel(5),
         axis.ticks.length = unit(5,"points") ,
         axis.title.y = element_text(size=15),
         axis.text = element_text(color = "black", size=15)
    )
  if (graph == T){
    ggsave(file = paste("BQ", paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}
data= read.table(file= "clipboard",sep= "\t",header =T)
mydata= read.table(file= "clipboard",sep= "\t",header =T)
mydata$genotype = factor(mydata$genotype, levels = c("Col-0", "npr1-1",  "npr4-4D", "n1n4")) # change based on your genotype names no spaces in genotype names
mydata$age = factor(mydata$age, levels = c("Y", "M"))
mydata$experiment = as.factor(mydata$experiment)
mydata$cfu = as.numeric(mydata$cfu)

YMBQGraph(mydata, barLabs = c("Col-0\n", "pen3-4\n", "pdr12-3\n", "pen3-4\npdr12-3"), ageCol = c( "#5D2169", "#FFFFFF"), exptID = "ARR-PEN3-all",
          #expCol = "#000000", 
          graph = F, box = F, width = 7, height = 6)

YMBQGraph(mydata[,], barLabs = c("Col-0\n", "npr1-1\n", "npr4-4D\n", "npr1-1\nnpr4-4D"), ageCol = c("#EEF3E2", "#00798C"), exptID = "ARR-NPR-23-2-only",expCol = "#000000", graph = F, box = F, width = 7, height = 6)
YMBQGraph(mydata[mydata$experiment == "ARR-NPR-23-2",], barLabs = c("Col-0\n", "npr1-1\n", "npr4-4D\n", "npr1-1\nnpr4-4D"), ageCol = c("#9609FF", "#FFFF00"), exptID = "ARR-NPR-23-2-only",expCol = "#000000", graph = F, box = F, width = 7, height = 6)
## deprecated
mydata = read.csv(file = file.choose(),header=T)
mydata= read.table(file= "clipboard",sep= "\t",header =T) ##This one copies a table from the clipboard
mydata = gather(mydata, genotype, measurement, Y_Col:M_b1b1)
mydata$measurement = log10(mydata$measurement)
mydata = mydata[!is.na(mydata$measurement), ]

mydata$genotype = factor(mydata$genotype, levels = c("Col-0", "npr1-1", "npr4-4D", "n1n4"))

mydata$genotype = factor(mydata$genotype, levels = c("Y_Col", "M_Col", "Y_ald1", "M_ald1",  "Y_sard4", "M_sard4", "Y_fmo1", "M_fmo1"))

mydata$genotype = factor(mydata$genotype, levels = c("Y_Col", "M_Col",  "Y_bak1.5", "M_bak1.5", "Y_bkk1.1", "M_bkk1.1", "Y_b1b1","M_b1b1"))

mydata$genotype = factor(mydata$genotype, levels = c("Y.Col", "M.Col",  "Y.pen3.4", "M.pen3.4", "Y.pdr12.3", "M.pdr12.3", "Y.p3p12","M.p3p12"))

mydata$genotype = factor(mydata$genotype, levels = c("Col-0", "bak1-5",  "bkk1-1", "b1b1"))
mydata$age = factor(mydata$age, levels = c("Y", "M"))
mydata$experiment = as.factor(mydata$experiment)
mydata$cfu = as.numeric(mydata$cfu)
# Plot
mydata %>%
  ggplot( aes(x=genotype, y=measurement )) +
  geom_boxplot(fill = rep(c("#595959","#FFFFFF"), 4)) + # the number of colours has to match the number of groups
  geom_jitter( width = 0.25, color=as.integer(as.factor(mydata$experiment)),
               size=2, alpha=0.5)  + coord_cartesian(ylim = c(4, 7.5))+scale_y_continuous(breaks=c(4, 5, 6, 7), 
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
