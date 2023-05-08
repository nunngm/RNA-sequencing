## Plotting qPCR data

library(tidyr)
library(ggplot2)
library(dplyr)
library(svglite)
library(agricolae)
library(car)
library(ggthemes)

setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-qRT-PCR\\Graphs")

df.avg = read.table("clipboard", sep = "\t", row.names = 1,header=T)

samps = rownames(df.avg)
samps = t(as.data.frame(strsplit(samps, split = "_", fixed = T)))
samps = as.data.frame(cbind(t(as.data.frame(strsplit(samps[,1], split = ".", fixed = T))), samps[,2]))
samps = samps[,1:2] #remove redundant info
colnames(samps) = c("genotype", "age")
samps$genotype = factor(samps$genotype, levels = c("Col","fmo1"))
samps$age = factor(samps$age, levels = c("Y", "M"))

df.graph = cbind(samps, df.avg)

upperbound <- function(x) {
  return(quantile(x, 0.75, na.rm = T) + 1.5 * IQR(x,na.rm = T))
}
lowerbound <- function(x){
  return(quantile(x, 0.25, na.rm = T) - 1.5 * IQR(x, na.rm=T))
}

is_outlier <- function(x) {
  return(x < quantile(x, 0.25, na.rm = T) - 1.5 * IQR(x, na.rm=T) | x > quantile(x, 0.75, na.rm = T) + 1.5 * IQR(x,na.rm = T))
}

targetGeneName = "RLP23"
refGeneName = "SEC5A"
data = df.graph %>%  mutate( sampGroup = paste(genotype, age, sep = "_"), .keep = "all") %>% mutate(target = 2^(-(get(targetGeneName)-get(refGeneName))))
anovaModel = aov(log2(target) ~ sampGroup, data = data)
print(HSD.test(anovaModel, alpha=0.05, "sampGroup", console=F)$groups)

## Graph for graphing 3 factor data of young and mature samples
qpcr3FGraph = function(data, targetGeneName, refGeneName, exptID = "temp", colours = c("red", "blue"), width = 8, height = 6, graph = F){
  data = data %>%  mutate( sampGroup = paste(genotype, age, sep = "_"), .keep = "all") %>% mutate(target = 2^(-(get(targetGeneName)-get(refGeneName))))
  print(data)
  anovaModel = aov(log2(target) ~ sampGroup, data = data)
  print(HSD.test(anovaModel, alpha=0.05, "sampGroup", console=F)$groups)
  
  p = data %>% group_by(sampGroup) %>% 
    #mutate(target = log10(target)) %>%
    #mutate(inlier = ifelse(is_outlier(target), as.numeric(NA), target), outlier = ifelse(is_outlier(target), target, as.numeric(NA)) ) %>% 
    mutate(inlier = target) %>% 
    ggplot(., aes(x=genotype, y=inlier, fill = age)) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 1), colour = "#000000", size = 0.75) +
    geom_jitter( size=2,#colour = df.graph$rep, 
                 alpha = 0.5, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.8)) +
    stat_summary(fun = mean,
                 fun.min = function(x) {ifelse(mean(x) - sd(x)>0, 
                                               mean(x) - sd(x)
                                               , 0 )
                   }, 
                 fun.max = function(x) {mean(x) + sd(x)}, 
                 geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000", position = position_dodge(width = 1)) +
    scale_y_continuous(expand = expansion(c(0, 0.1)))+
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.x =  element_blank()
    ) + ylab(paste0(targetGeneName,"/", refGeneName)) + xlab("") +
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          #strip.text.x = element_text(size = 15),
          #strip.background = element_rect(colour = "black", fill = "#FFFFFF", size = 1),
          axis.line = element_line(colour = "black", size=0),
          axis.title.x=element_text(size=15),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length = unit(5,"points") ,
          axis.title.y = element_text(size=15),
          axis.text = element_text(color = "black", size=15),
          strip.background.x = element_blank(),
          strip.text.x = element_blank()) +
    scale_fill_manual(values = colours)
  # + theme_few()
  if(graph ==T){
    exptID = readline(prompt = "Enter experimentID:")
    ggsave(file = paste(targetGeneName, refGeneName, paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}

qpcr3FGraph(df.graph, targetGeneName = "RLP28", refGeneName = "CUL4",exptID = "ARR-PIP-22-1-fmo1", height = 6, width = 7, colours = c("#F6A63C", "#54B031"), graph = T)


targetGeneName = "RLP28"
refGeneName = "SEC5A"
qpcr3FGraph(df.graph, targetGeneName = targetGeneName, refGeneName = refGeneName,exptID = "ARR-PIP-22-1", height = 6, width = 7, colours = c("#54B031", "#0993AE" , "#F6A63C"), graph = F)
data = df.graph %>%  mutate( sampGroup = paste(age, treatment, hpi, sep = "_"), .keep = "all") %>% mutate(target = 2^(-(get(targetGeneName)-get(refGeneName))))
anovaModel = aov(log2(target) ~ sampGroup, data = data[data$age=="M",])
HSD.test(anovaModel, alpha=0.05, "sampGroup", console=F)$groups
anovaModel = aov(log2(target) ~ sampGroup, data = data)
temp = HSD.test(anovaModel, alpha=0.05, "sampGroup", console=F)$groups
temp[,1] = 2^temp[,1]
temp

qpcrCtGraph = function(data, targetGeneName, exptID = "exptID", colours = c("red", "green", "blue"), width = 8, height = 6, graph = F){
  data = data %>%  mutate(target = get(targetGeneName), sampGroup = paste(genotype, age, sep = "_"), .keep = "all")
  print(data)
  anovaModel = aov(target ~ sampGroup, data = data)
  print(HSD.test(anovaModel, alpha=0.05, "sampGroup", console=F)$groups)
  
  p = data %>% group_by(sampGroup) %>% 
    mutate(inlier = ifelse(is_outlier(target), as.numeric(NA), target), outlier = ifelse(is_outlier(target), target, as.numeric(NA)) ) %>%
    ggplot(., aes(x=genotype, y=inlier, fill = age)) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 1), colour = "#000000", size = 0.75) +
    geom_jitter( size=2, alpha = 0.5, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.8)) + 
    stat_summary(fun = mean,
                 fun.min = function(x) {ifelse(mean(x) - sd(x)>0,mean(x) - sd(x),0 )}, 
                 fun.max = function(x) {mean(x) + sd(x)}, 
                 geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000", position = position_dodge(width = 1)) +
    scale_y_continuous(expand = expansion(c(0, 0.1)))+
    theme(
      legend.position="right",
      plot.title = element_text(size=11),
      axis.text.x =  element_blank()
    ) + ylab(paste0(targetGeneName,"(Ct)")) + xlab("") +
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          strip.text.x = element_text(size = 15),
          strip.background = element_rect(colour = "black", fill = "#FFFFFF", size = 1),
          axis.line = element_line(colour = "black", size=0),
          axis.title.x=element_text(size=15),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length = unit(5,"points") ,
          axis.title.y = element_text(size=15),
          axis.text = element_text(color = "black", size=15)) +
    scale_fill_manual(values = colours)
  # + theme_few()
  if(graph ==T){
    ggsave(file = paste("CtVal", targetGeneName, paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}

qpcrCtGraph(df.graph, targetGeneName = "CUL4", exptID = "ARR-PIP-22-1-genotype", colours = c("#54B031", "#0993AE" ), graph = T)


# Moving on to weekly expression graphs

weeklyBacterialLevel = function(data, exptID = "exptID", colours = c("red", "green", "blue"), width = 8, height = 6, graph = F){
  data = data %>%  mutate(target = 2^(-(get(targetGeneName)-get(refGeneName))), sampGroup = paste(age, treatment, hpi, sep = "_"), .keep = "all")
  print(data)
  anovaModel = aov(log2(target) ~ sampGroup, data = data)
  print(HSD.test(anovaModel, alpha=0.05, "sampGroup", console=F)$groups)
  
  p = data %>% group_by(sampGroup) %>% 
    mutate(inlier = ifelse(is_outlier(target), as.numeric(NA), target), outlier = ifelse(is_outlier(target), target, as.numeric(NA)) ) %>%
    ggplot(., aes(x=hpi:treatment, y=inlier, fill = treatment)) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 1), colour = "#000000", size = 0.75) +
    geom_jitter( size=2, alpha = 0.5, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.8)) + facet_grid(.~age, labeller = labeller(age = c(Y = "Young", M = "Mature"))) +
    stat_summary(fun = mean,
                 fun.min = function(x) {ifelse(mean(x) - sd(x)>0,mean(x) - sd(x),0 )}, 
                 fun.max = function(x) {mean(x) + sd(x)}, 
                 geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000", position = position_dodge(width = 1)) +
    scale_y_continuous(expand = expansion(c(0, 0.1)))+
    theme(
      legend.position="right",
      plot.title = element_text(size=11),
      axis.text.x =  element_blank()
    ) + ylab(paste0(targetGeneName,"/", refGeneName)) + xlab("") +
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          strip.text.x = element_text(size = 15),
          strip.background = element_rect(colour = "black", fill = "#FFFFFF", size = 1),
          axis.line = element_line(colour = "black", size=0),
          axis.title.x=element_text(size=15),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length = unit(5,"points") ,
          axis.title.y = element_text(size=15),
          axis.text = element_text(color = "black", size=15)) +
    scale_fill_manual(values = colours)
  # + theme_few()
  if(graph ==T){
    ggsave(file = paste(targetGeneName, refGeneName, paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}



qpcrWeeklyGraph = function(data, targetGeneName, refGeneName, exptID = "exptID", colours = c("red", "green", "blue"), width = 8, height = 6, graph = F){
  data = data %>%  mutate(target = 2^(-(get(targetGeneName)-get(refGeneName))), sampGroup = paste(age, treatment, hpi, sep = "_"), .keep = "all")
  print(data)
  anovaModel = aov(log2(target) ~ sampGroup, data = data)
  print(HSD.test(anovaModel, alpha=0.05, "sampGroup", console=F)$groups)
  
  p = data %>% group_by(sampGroup) %>% 
    mutate(inlier = ifelse(is_outlier(target), as.numeric(NA), target), outlier = ifelse(is_outlier(target), target, as.numeric(NA)) ) %>%
    ggplot(., aes(x=hpi:treatment, y=inlier, fill = treatment)) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 1), colour = "#000000", size = 0.75) +
    geom_jitter( size=2, alpha = 0.5, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.8)) + facet_grid(.~age, labeller = labeller(age = c(Y = "Young", M = "Mature"))) +
    stat_summary(fun = mean,
                 fun.min = function(x) {ifelse(mean(x) - sd(x)>0,mean(x) - sd(x),0 )}, 
                 fun.max = function(x) {mean(x) + sd(x)}, 
                 geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000", position = position_dodge(width = 1)) +
    scale_y_continuous(expand = expansion(c(0, 0.1)))+
    theme(
      legend.position="right",
      plot.title = element_text(size=11),
      axis.text.x =  element_blank()
    ) + ylab(paste0(targetGeneName,"/", refGeneName)) + xlab("") +
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          strip.text.x = element_text(size = 15),
          strip.background = element_rect(colour = "black", fill = "#FFFFFF", size = 1),
          axis.line = element_line(colour = "black", size=0),
          axis.title.x=element_text(size=15),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length = unit(5,"points") ,
          axis.title.y = element_text(size=15),
          axis.text = element_text(color = "black", size=15)) +
    scale_fill_manual(values = colours)
  # + theme_few()
  if(graph ==T){
    ggsave(file = paste(targetGeneName, refGeneName, paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}

p = df.graph %>%
  group_by(age, hpi, treatment) %>% 
  mutate(target = 2^(-(get(targetGeneName)-get(refGeneName))) , .keep = "unused") %>% 
  mutate(inlier = ifelse(is_outlier(target), as.numeric(NA), target), outlier = ifelse(is_outlier(target), target, as.numeric(NA)) ) %>%
ggplot(., aes(x=hpi:treatment, y=inlier, fill = treatment)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 1), colour = "#000000", size = 0.75) +
  stat_summary(fun = mean,
               fun.min = function(x) {ifelse(mean(x) - sd(x)>0,mean(x) - sd(x),0 )}, 
               fun.max = function(x) {mean(x) + sd(x)}, 
               geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000", position = position_dodge(width = 1)) +
  geom_jitter( size=2, alpha=1, position = position_jitterdodge(dodge.width = 1, jitter.width = 0)) + facet_grid(.~age) +
  scale_y_continuous(expand = expansion(c(0, 0.1)))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + ylab(paste0(targetGeneName,"/SEC5A")) + theme_few()
  geom_text(aes(label=Tukey, y = uptake_mean + sd + 2), size = 3, color = "Gray25",
            show.legend = FALSE,
            position = position_dodge(0.9))


p = expression %>%
  group_by(type) %>%
  mutate(inlier = ifelse(is_outlier(!!as.name(targetGeneName)), as.numeric(NA), !!as.name(targetGeneName)), outlier = ifelse(is_outlier(!!as.name(targetGeneName)), !!as.name(targetGeneName), as.numeric(NA)) ) %>%
  ggplot(., aes(x=type, y=inlier, colour = rep)) +
  stat_summary(fun = mean, geom = "bar", fill = rep(c( "#444444",
                                                       "#666666", "#9A9A9A", "#CDCDCD", "#FFFFFF"),2), colour = "#000000", size = 0.75) +
  stat_summary(fun = mean,
               fun.min = function(x) {mean(x) - sd(x)}, 
               fun.max = function(x) {mean(x) + sd(x)}, 
               geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000") +
  #geom_boxplot(fill = rep(c("#FFFFFF"), 5)) +
  geom_jitter(width = 0.25, color= "#000000", size = 2, alpha = 0.4) +
  geom_point(aes(x = type, y = outlier), size =2, alpha = 1, shape = 8, colour = "#000000") +
  scale_y_continuous(expand = expansion(c(0, 0.1)))   +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  xlab("Weeks post-germination (wpg)") + ylab(paste0(targetGeneName,"/SEC5A")) + theme(panel.grid.major = element_blank(),  
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
p 
ggsave(file = paste0(targetGeneName,"_WKEX-22-2.svg"), plot = p, width = 5, height = 4)