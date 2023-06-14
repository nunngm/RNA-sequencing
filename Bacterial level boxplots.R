## Making Boxplots

# Libraries
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(grDevices)
library(agricolae)

#laptop directory
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")

#Desktop directory
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")

mydata = read.csv(file = file.choose(),header=T)
mydata= read.table(file= "clipboard",sep= "\t",header =T)
mydata = gather(mydata, genotype, measurement, Y_Col:M_b1b1)
mydata$measurement = log10(mydata$measurement)
mydata = mydata[!is.na(mydata$measurement), ]

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

YMBQGraph = function(data, ageCol = c("#595959","#FFFFFF"), ylim = c(4,8), expCol = NA, graph = F, width = 5, height = 4, exptID = "temp", box = F){
  data = data %>% mutate(sampGroup = paste(age, genotype, sep = "_"), cfu = log10(cfu), .keep = "all")
  if(is.na(expCol)){expCol = as.integer(as.factor(data$experiment))}
  
  faces = c("plain", rep("italic", times = length(levels(data$genotype))-1))
  print(data)
  anovaModel = aov(cfu ~ sampGroup, data = data)
  print(data %>%group_by(genotype,age)%>% summarise(cfu = mean(cfu)) %>% summarise(foldDiff = 10^(cfu[age == 'Y'] - cfu[age == 'M'])))
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
                       expand = c(0, 0.05))   + #scale_x_discrete(labels = barLabs) +
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

mydata= read.table(file= "clipboard",sep= "\t",header =T)
mydata$genotype = factor(mydata$genotype, levels = c("Col-0", "ald1-T2",  "fmo1-1"))
mydata$age = factor(mydata$age, levels = c("Y", "M"))
mydata$experiment = as.factor(mydata$experiment)
mydata$cfu = as.numeric(mydata$cfu)

YMBQGraph(temp, exptID = "ARR-CSR-23-1",expCol = "#000000", graph = T, box = F, width = 7, height = 5)

