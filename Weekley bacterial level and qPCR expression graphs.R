## Weekly expression and bacterial level plots

library(tidyr)
library(ggplot2)
library(dplyr)
library(svglite)
library(agricolae)
library(car)
library(ggthemes)
library(grDevices)

#set working directory and read in the data
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-wkly bac qRT-PCR\\output")
df.avg = read.table("clipboard", sep = "\t", row.names = 1, header=T)

upperbound <- function(x) {
  return(quantile(x, 0.75, na.rm = T) + 1.5 * IQR(x,na.rm = T))
}
lowerbound <- function(x){
  return(quantile(x, 0.25, na.rm = T) - 1.5 * IQR(x, na.rm=T))
}

is_outlier <- function(x) {
  return(x < quantile(x, 0.25, na.rm = T) - 1.5 * IQR(x, na.rm=T) | x > quantile(x, 0.75, na.rm = T) + 1.5 * IQR(x,na.rm = T))
}
# collect and organise sample information
samps = rownames(df.avg)
samps = as.data.frame(t(as.data.frame(strsplit(samps, split = "_", fixed = T))))
#samps = as.data.frame(cbind(t(as.data.frame(strsplit(samps[,1], split = ".", fixed = T))), samps[,2]))
colnames(samps) = c("age","rep")
samps$age = factor(samps$age)

df.graph = cbind(samps, df.avg)
pal = colorRampPalette(c( "white", "#54B031"))

weeklyBacterialLevel = function(data, labels = NULL, exptID = "exptID", colours = c("yellow","red", "green", "blue"), width = 5, height = 4, graph = F){
  data = data %>%  mutate(cfu = log10(cfu), .keep = "all")
  print(data)
  anovaModel = aov(cfu ~ age, data = data)
  print(HSD.test(anovaModel, alpha=0.05, "age", console=F)$groups)
  
  p = data %>% group_by(age) %>% 
    ggplot(., aes(x=age, y=cfu, fill = age)) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 1), colour = "#000000",width =1, size = 1) +
    geom_jitter( size=2, alpha = 0.5, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.8)) +
    stat_summary(fun = mean,
                 fun.min = function(x) {ifelse(mean(x) - sd(x)>0,mean(x) - sd(x),0 )}, 
                 fun.max = function(x) {mean(x) + sd(x)}, 
                 geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000", position = position_dodge(width = 1)) +
    scale_y_continuous(breaks=c(4, 5, 6, 7), 
                       labels = expression(10^4, 10^5, 10^6, 10^7),
                       expand = c(0, 0)) +coord_cartesian(ylim = c(4, 7.5)) +
    scale_x_discrete( labels = labels) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +ylab(bquote('Bacterial level (cfu leaf disc'^-1*')')) + xlab("") +
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_blank(),
          strip.text.x = element_text(size = 15),
          strip.background = element_rect(colour = "black", fill = "#FFFFFF", size = 1),
          axis.line = element_line(colour = "black", size=1),
          axis.title.x=element_text(size=15),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length = unit(5,"points") ,
          axis.title.y = element_text(size=15),
          axis.text = element_text(color = "black", size=20)) +
    scale_fill_manual(values = colours)
  
  # + theme_few()
  if(graph ==T){
    ggsave(file = paste(exptID, "cfu.svg", sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}

weeklyBacterialLevel(data = df.graph, labels = c(3,4,5,6), exptID = "WKEX-22-3", colours = pal(4), graph = T)

## weekly qPCR data
weeklyQPCR = function(data, targetGeneName, scale = c(0,NA), refGeneName, labels = NULL, exptID = "exptID", colours = c("yellow","red", "green", "blue"), width = 5, height = 4, graph = F){
  data = data %>%  mutate(target = 2^(-(get(targetGeneName)-get(refGeneName))), .keep = "all")
  print(data)
  anovaModel = aov(log2(target) ~ age, data = data)
  print(HSD.test(anovaModel, alpha=0.05, "age", console=F)$groups)
  
  p = data %>% group_by(age) %>% mutate(inlier = ifelse(is_outlier(target), as.numeric(NA), target), outlier = ifelse(is_outlier(target), target, as.numeric(NA)) ) %>%
    ggplot(., aes(x=age, y=inlier, fill = age)) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 1), colour = "#000000",width =1, size = 1) +
    geom_jitter( size=2, alpha = 0.5, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.8)) +
    stat_summary(fun = mean,
                 fun.min = function(x) {ifelse(mean(x) - sd(x)>0,mean(x) - sd(x),0 )}, 
                 fun.max = function(x) {mean(x) + sd(x)}, 
                 geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000", position = position_dodge(width = 1)) +
       geom_point(aes(x = age, y = outlier), size =2, alpha = 1, shape = 17, colour = "#000000") +
    scale_y_continuous(limits = scale, expand = expansion(c(0, 0.1))) +
    scale_x_discrete( labels = labels) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +ylab(paste0(targetGeneName,"/", refGeneName)) + xlab("") +
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_blank(),
          strip.text.x = element_text(size = 15),
          strip.background = element_rect(colour = "black", fill = "#FFFFFF", size = 1),
          axis.line = element_line(colour = "black", size=1),
          axis.title.x=element_text(size=15),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length = unit(5,"points") ,
          axis.title.y = element_text(size=15),
          axis.text = element_text(color = "black", size=20)) +
    scale_fill_manual(values = colours)
  
  # + theme_few()
  if(graph ==T){
    ggsave(file = paste(exptID, targetGeneName, refGeneName, "qPCR.svg", sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}

weeklyQPCR(data = df.graph, scale = c(0,0.16), targetGeneName = "FMO1", refGeneName = "SEC5A", labels = c(3,4,5,6), colours = pal(4), exptID = "WKEX-23-2-scaled" , graph = T)
# outlier calculation
