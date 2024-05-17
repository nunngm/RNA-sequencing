##Quantification of SA levels
library(ggplot2)
library(dplyr)
library(agricolae)
library(tidyr)
library(stats)

## From raw plate
layout = unlist(read.table(file = "clipboard", sep = "\t", header = F))
lum = unlist(read.table(file = "clipboard", sep = "\t", header = F))

names(layout) == names(lum)
linData = as.data.frame(cbind(names(layout), layout, lum))

curveData = linData[grepl("^{0,1}[0-9]{0,}.{0,1}[0-9]{1,}$", linData$layout),]
colnames(curveData) = c("ID", "SA", "lum")
curveData$SA = as.numeric(curveData$SA)
curveData$lum = as.numeric(curveData$lum)

## Fitting to a model
# plot(curveData$layout, curveData$lum)
# fit = nls(lum ~ SSfpl(SA, A, B, xmid, scal), data = curveData )
# fit = nls( lum ~ SSlogis(SA, Asym, xmid, scal), curveData)
# fit = nls(lum ~ SSgompertz(SA, Asym, b2, b3), data = curveData)

### an exponential model worked best
 # curveData = curveData[curveData$SA != 0,] # For a log transformation there can be no 0 values
 # fit.exp = lm(log(SA) ~ lum, data = curveData) # the model

### No actually a linear model was by far the best.
curveData = curveData[curveData$SA != 200,]
fit = lm(SA ~ lum, data = curveData) #kill me, I spent so long learning non-linear fits just for a linear fit to be


# summary(fit)
# summary(fit.exp) # check the model parameters
summary(fit)


# new.data = data.frame(SA = seq(0, 200, length.out = 24))
# predict(fit, curveData$lum)
# interval = as_tibble(predFit(fit, newdata = new.data, interval = "confidence", level = 0.99)) %>% mutate(SA = new.data$SA)
# interval$value = exp(interval$value)
# p1 <- ggplot(curveData) +  geom_point(aes(x=SA, y=lum),size=2, colour="black") + xlab("Time (h)") + ylab("Optical density (OD600)") 
# p1 + geom_line(data=interval, aes(x = SA, y = fit )) #+
#   geom_ribbon(data=interval, aes(x=SA, ymin=lwr, ymax=upr), alpha=0.5, inherit.aes=F, fill="blue")+
#   theme_classic()

new.data = data.frame(lum = seq(0, max(curveData$lum), length.out = 22))
interval = as_tibble(predFit(fit, newdata = new.data)) %>% mutate(lum = new.data$lum)
# interval$value = exp(interval$value)
p1 <- ggplot(curveData) +  geom_point(aes(x=lum, y=SA),size=2, colour="black") + xlab("Luminescence") + ylab("ng SA") 
p1 + geom_line(data=interval, aes(x = lum, y = value )) +  theme_classic()

#linData$predictedSA = exp(predict(fit.exp, newdata =  data.frame(lum = as.numeric(linData$lum))))
linData$predictedSA = predict(fit, newdata =  data.frame(lum = as.numeric(linData$lum)))
linData = linData[linData$layout != "EMPTY" & !grepl("^{0,1}[0-9]{0,}.{0,1}[0-9]{1,}$", linData$layout),]

linData = cbind(linData, t(as.data.frame(strsplit(linData$layout, split = "_"), row.names =  c("cellType", "hpi"))))



write.table(linData, "clipboard", sep = "\t", row.names = F)
        
# From summarized data
mydata= read.table(file= "clipboard",sep= "\t",header =T)
mydata$genotype = factor(mydata$genotype, levels = c("Col-0", "pen3-4", "pdr12-3", "p3p12"))
mydata$experiment = as.factor(mydata$experiment)
mydata$hpi = factor(mydata$hpi, levels = c("24"))
mydata$treatment = factor(mydata$treatment, levels = c("UN", "MO", "PST"))
mydata$intercellular = as.numeric(mydata$intercellular)
mydata$intracellular = as.numeric(mydata$intracellular)


SAQuantGraph = function(data, genotypeCol = c("#378717","#FFFF00", "#6DFDFD"), ylim = c(4,8), expCol = NA, graph = F, width = 5, height = 4, exptID = "temp", box = F){
}
SAQuantGraph = function(data, genotypes = c("Col-0", "ald1-T2", "fmo1-1"), colours = c("#378717","#FFFF00", "#6DFDFD"), selectTimes = c("UN", "12", "18", "24"), 
                        barLabs = c("Untreated", "12 hpi", "18 hpi", "24 hpi"), ylim = c(0,NA), width = 5, height = 4, SA = c("inter", "intra"),
                        graph = F){
  if (length(SA) == 2){
    SA = "inter"
  }
  if (SA == "inter"){
    selectColumn = "intercellular"
    lab.y = bquote('Intercellular SA (ng ml IWF'^-1*')')
  } else{
    selectColumn = "intracellular"
    lab.y = bquote('Intracellular SA (ng gfw'^-1*')')
  }
  df = data[data$hpi %in% selectTimes, ] #Rows at the selected times
  df$genotype = factor(df$genotype, levels = genotypes)
  df$experiment = as.factor(df$experiment)
  df$hpi = factor(df$hpi, levels = selectTimes)
  df$hpi = factor(df$hpi, levels = selectTimes)
  df$intercellular = as.numeric(df$intercellular)
  df$intracellular = as.numeric(df$intracellular)
  #take full data and pare down to just the required treatment and untreated
  
  p = ggplot(df, aes(x = hpi, y = get(selectColumn), fill = genotype, group = genotype)) + 
    stat_summary(aes(group = genotype), colour = "#000000", fun = mean, geom = 'bar', width = 0.6, size = 1, position = position_dodge(width = 0.8)) + 
    stat_summary( aes(y = get(selectColumn), group = genotype), fun = mean,
                  fun.min = function(x) {pmax(mean(x) - sd(x), 0, na.rm = T)}, 
                  fun.max = function(x) {mean(x) + sd(x)}, 
                  geom = "errorbar", lty =1 , size =0.75, width = 0.25, position = position_dodge(width = 0.8)) +
    geom_jitter( alpha = 0.4, size = 2, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8)) +
    xlab("") +   scale_y_continuous(limits = ylim, expand = c(0,0)) + ylab(lab.y) +  scale_x_discrete(labels = barLabs) +
    scale_fill_manual(values = colours) +
    #ylim(NA,3.2)  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black", size=1),
          axis.title.x=element_blank(),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length.x = unit(5, "points"),
          axis.ticks.length.y = unit(5,"points") ,
          axis.title.y = element_blank(),
          axis.text = element_text(colour = "black", size=15),
          axis.text.x = element_text(vjust = -0.25),
          legend.position = "none",
          plot.margin = unit(c(20,0,10,0), "points")
    )
  #Full anova
  data$genohpi = paste(data$genotype,data$hpi, sep = "_")
  anovaModel = aov(get(selectColumn) ~ genohpi, data = data)
  print(HSD.test(anovaModel, alpha=0.05, "genohpi", console=F)$groups)
  
  if(graph ==T){
    exptID = readline(prompt = "Enter the experiment ID: ")
    ggsave(file = paste(selectColumn, "SAQuant", paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}

SAQuantGraph.byTreatment = function(data, genotypes = c("Col-0", "pen3-4", "pdr12-3", "p3p12"), colours = c("#FFFFFF", "#00BBFF", "#FF3853", "#5D2169"), selectTreatments = c("UN", "MO", "PST"),
                        barLabs = c("Col-0\n", "pen3-4\n", "pdr12-3\n", "pen3-4\npdr12-3"), ylim = c(0,NA), width = 5, height = 4, SA = c("inter", "intra"),
                        graph = F){
  if (length(SA) == 2){
    SA = "inter"
  }
  if (SA == "inter"){
    selectColumn = "intercellular"
    lab.y = bquote('Intercellular SA (ng ml IWF'^-1*')')
  } else{
    selectColumn = "intracellular"
    lab.y = bquote('Intracellular SA (ng gfw'^-1*')')
  }
  df = data[data$treatment %in% selectTreatments, ] #Rows at the selected times
  df$genotype = factor(df$genotype, levels = genotypes)
  df$experiment = as.factor(df$experiment)
  df$treatment = factor(df$treatment, levels = selectTreatments)
  df$intercellular = as.numeric(df$intercellular)
  df$intracellular = as.numeric(df$intracellular)
  faces = c("plain", rep("italic", times = length(levels(data$genotype))-1))
  #take full data and pare down to just the required treatment and untreated
  
  p = ggplot(df, aes(x = genotype, y = get(selectColumn), fill = treatment, group = treatment)) + 
    stat_summary(aes(group = treatment), colour = "#000000", fun = mean, geom = 'bar', width = 0.6, size = 1, position = position_dodge(width = 0.8)) + 
    stat_summary( aes(y = get(selectColumn), group = treatment), fun = mean,
                  fun.min = function(x) {pmax(mean(x) - sd(x), 0, na.rm = T)}, 
                  fun.max = function(x) {mean(x) + sd(x)}, 
                  geom = "errorbar", lty =1 , size =0.75, width = 0.25, position = position_dodge(width = 0.8)) +
    geom_jitter( alpha = 0.4, size = 2, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8)) +
    xlab("") +   scale_y_continuous(limits = ylim, expand = c(0,0)) + ylab(lab.y) +  scale_x_discrete(labels = barLabs) +
    scale_fill_manual(values = colours) +
    #ylim(NA,3.2)  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black", size=1),
          axis.title.x=element_blank(),
          #axis.text.x=element_blank()),
          axis.text.x=element_text(face = faces),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length.x = unit(5, "points"),
          axis.ticks.length.y = unit(5,"points") ,
          axis.title.y = element_blank(),
          axis.text = element_text(colour = "black", size=20),
          legend.position = "none",
          plot.margin = unit(c(20,0,10,0), "points")
    )
  #Full anova
  data$genotreat = paste(data$genotype,data$treatment, sep = "_")
  anovaModel = aov(get(selectColumn) ~ genotreat, data = data)
  print(HSD.test(anovaModel, alpha=0.05, "genotreat", console=F)$groups)
  
  if(graph ==T){
    exptID = readline(prompt = "Enter the experiment ID: ")
    ggsave(file = paste(selectColumn, "SAQuant", paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}
SAQuantGraph.byTreatment(mydata[mydata$age == "M" & mydata$experiment =="ARR-PEN3-21-2",], ylim = c(0,4000), selectTreatments = c("UN","PST"), colours = c("#FFFFFF", "#FF3853"), width =7, height = 6, SA="intra", graph = T)
## Working on Pst strains


mydata= read.table(file= "clipboard",sep= "\t",header =T)
df = mydata
df$genotype = factor(df$genotype, levels = c("untreated","mock","Pst", "Pst ΔalgU mucAB",  "Pst ΔalgU mucAB ΔalgD"))

df$experiment = as.factor(df$experiment)
df$hpi = factor(df$hpi, levels = c("6", "12", "24"))
df$intercellular = as.numeric(df$intercellular)
df$genohpi = paste(df$genotype,df$hpi, sep = "_")
df$genohpi = factor(df$genohpi, levels = c("Pst_6","Pst_12", "Pst_24", "Pst ΔalgU mucAB_6", "Pst ΔalgU mucAB_12", "Pst ΔalgU mucAB_24", "Pst ΔalgU mucAB ΔalgD_6", "Pst ΔalgU mucAB ΔalgD_12", "Pst ΔalgU mucAB ΔalgD_24"))

df$genohpi = paste(df$genotype,df$treatment,df$hpi, df$post, sep = "_")
df$genohpi = factor(df$genohpi, levels = c("Col-0_mock_6_hpt", "Col-0_mock_12_hpt", "Col-0_mock_24_hpt", "Col-0_mock_6_hpi", "Col-0_mock_12_hpi", "Col-0_mock_24_hpi", "Col-0_flg22_6_hpt", "Col-0_flg22_12_hpt", "Col-0_flg22_24_hpt", "Col-0_flg22_6_hpi", "Col-0_flg22_12_hpi", "Col-0_flg22_24_hpi",
                                           "fls2_mock_6_hpt", "fls2_mock_12_hpt", "fls2_mock_24_hpt", "fls2_mock_6_hpi", "fls2_mock_12_hpi", "fls2_mock_24_hpi", "fls2_flg22_6_hpt", "fls2_flg22_12_hpt", "fls2_flg22_24_hpt", "fls2_flg22_6_hpi", "fls2_flg22_12_hpi", "fls2_flg22_24_hpi"
                                           ))

p = ggplot(df, aes(x = genohpi, y = intercellular, fill= "white")) + 
  stat_summary( colour = "#000000", fun = mean, geom = 'bar', width = 0.6, size = 1, position = position_dodge(width = 0.8))  +
  stat_summary( fun = mean,
                fun.min = function(x) {pmax(mean(x) - sd(x), 0, na.rm = T)}, 
                fun.max = function(x) {mean(x) + sd(x)}, 
                geom = "errorbar", lty =1 , size =0.75, width = 0.25, position = position_dodge(width = 0.8))  +
  #geom_jitter( alpha = 0.4, size = 2 ,width = 0, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8))+
  #) +
  xlab("") + scale_y_continuous(limits = c(0,2500), expand = c(0,0)) + ylab("") +  scale_x_discrete(labels = rep(c("6", "12", "24"), times = 8))  +
  scale_fill_manual(values = "white") +
  #ylim(NA,3.2)  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size=1),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank()),
        axis.ticks=element_line(colour = "black", size =1),
        axis.ticks.length.x = unit(5, "points"),
        axis.ticks.length.y = unit(5,"points") ,
        axis.title.y = element_text(size=15),
        axis.text = element_text(colour = "black", size=18),
        axis.text.x = element_text(vjust = -0.25),
        legend.position = "none",
        plot.margin = unit(c(20,0,10,0), "points")
  )
p
mydata$combined = paste(mydata$age, mydata$genotype, mydata$hpi, sep = "_")
anovaModel = aov(log2(intracellular+1) ~ combined, data = temp)
print(HSD.test(anovaModel, alpha=0.05, "combined", console=F)$groups)

ggsave(file = "PTI-IWF-20-3.svg", plot = p, width = 16, height = 8)


setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-PTI&Biofilms paper\\Exp-Triple and Quad\\IWF SA levels")
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR-NPR&PEN3\\exp-PEN3 SA levels\\graphs")
SAQuantGraph(mydata, genotypes = c("Pst", "Pst ΔalgU mucAB",  "Pst ΔalgU mucAB ΔalgD"), selectTimes = c("6", "12", "24"), barLabs = c("6 hpi", "12 hpi", "24 hpi"), SA = "inter")
SAQuantGraph(mydata, genotypes = c("Col-0", "pen3-4", "pdr12-3", "p3p12"), colours = c("#FFFFFF", "#00BBFF", "#FF3853", "#5D2169"), selectTimes = c("UN","24"), barLabs = c("Untreated"), SA = "inter", ylim = c(0,2500), graph = T, height = 6, width = 8)
SAQuantGraph.byTreatment(mydata, genotypes)
#Double mutant colours
colours = c("#FFFFFF", "#00BBFF", "#FF3853", "#5D2169")

# Working with yeast quantification
SAQuantGraph = function(data, genotypes = c("Col-0", "ald1-T2", "fmo1-1"), colours = c("#378717","#FFFF00", "#6DFDFD"), selectTimes = c("UN", "12", "18", "24"), 
                        barLabs = c("Untreated", "12 hpi", "18 hpi", "24 hpi"), ylim = c(0,NA), width = 5, height = 4, SA = c("inter", "intra"),
                        graph = F){
  if (length(SA) == 2){
    SA = "inter"
  }
  if (SA == "inter"){
    selectColumn = "intercellular"
    lab.y = bquote('Intercellular SA (ng ml IWF'^-1*')')
  } else{
    selectColumn = "intracellular"
    lab.y = bquote('Intracellular SA (ng gfw'^-1*')')
  }
  df = data[data$hpi %in% selectTimes, ] #Rows at the selected times
  df$genotype = factor(df$genotype, levels = genotypes)
  df$experiment = as.factor(df$experiment)
  df$hpi = factor(df$hpi, levels = selectTimes)
  df$hpi = factor(df$hpi, levels = selectTimes)
  df$intercellular = as.numeric(df$intercellular)
  df$intracellular = as.numeric(df$intracellular)
  #take full data and pare down to just the required treatment and untreated
  
  p = ggplot(df, aes(x = hpi, y = get(selectColumn), fill = genotype, group = genotype)) + 
    stat_summary(aes(group = genotype), colour = "#000000", fun = mean, geom = 'bar', width = 0.6, size = 1, position = position_dodge(width = 0.8)) + 
    stat_summary( aes(y = get(selectColumn), group = genotype), fun = mean,
                  fun.min = function(x) {pmax(mean(x) - sd(x), 0, na.rm = T)}, 
                  fun.max = function(x) {mean(x) + sd(x)}, 
                  geom = "errorbar", lty =1 , size =0.75, width = 0.25, position = position_dodge(width = 0.8)) +
    geom_jitter( alpha = 0.4, size = 2, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8)) +
    xlab("") +   scale_y_continuous(limits = ylim, expand = c(0,0)) + ylab(lab.y) +  scale_x_discrete(labels = barLabs) +
    scale_fill_manual(values = colours) +
    #ylim(NA,3.2)  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black", size=1),
          axis.title.x=element_blank(),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length.x = unit(5, "points"),
          axis.ticks.length.y = unit(5,"points") ,
          axis.title.y = element_blank(),
          axis.text = element_text(colour = "black", size=15),
          axis.text.x = element_text(vjust = -0.25),
          legend.position = "none",
          plot.margin = unit(c(20,0,10,0), "points")
    )
  #Full anova
  data$genohpi = paste(data$genotype,data$hpi, sep = "_")
  anovaModel = aov(get(selectColumn) ~ genohpi, data = data)
  print(HSD.test(anovaModel, alpha=0.05, "genohpi", console=F)$groups)
  
  if(graph ==T){
    exptID = readline(prompt = "Enter the experiment ID: ")
    ggsave(file = paste(selectColumn, "SAQuant", paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}

SAQuantGraph.byTreatment = function(data, genotypes = c("Col-0", "pen3-4", "pdr12-3", "p3p12"), colours = c("#FFFFFF", "#00BBFF", "#FF3853", "#5D2169"), selectTreatments = c("UN", "MO", "PST"),
                                    barLabs = c("Col-0\n", "pen3-4\n", "pdr12-3\n", "pen3-4\npdr12-3"), ylim = c(0,NA), width = 5, height = 4, SA = c("inter", "intra"),
                                    graph = F){
  if (length(SA) == 2){
    SA = "inter"
  }
  if (SA == "inter"){
    selectColumn = "intercellular"
    lab.y = bquote('Intercellular SA (ng ml IWF'^-1*')')
  } else{
    selectColumn = "intracellular"
    lab.y = bquote('Intracellular SA (ng gfw'^-1*')')
  }
  df = data[data$treatment %in% selectTreatments, ] #Rows at the selected times
  df$genotype = factor(df$genotype, levels = genotypes)
  df$experiment = as.factor(df$experiment)
  df$treatment = factor(df$treatment, levels = selectTreatments)
  df$intercellular = as.numeric(df$intercellular)
  df$intracellular = as.numeric(df$intracellular)
  faces = c("plain", rep("italic", times = length(levels(data$genotype))-1))
  #take full data and pare down to just the required treatment and untreated
  
  p = ggplot(df, aes(x = genotype, y = get(selectColumn), fill = treatment, group = treatment)) + 
    stat_summary(aes(group = treatment), colour = "#000000", fun = mean, geom = 'bar', width = 0.6, size = 1, position = position_dodge(width = 0.8)) + 
    stat_summary( aes(y = get(selectColumn), group = treatment), fun = mean,
                  fun.min = function(x) {pmax(mean(x) - sd(x), 0, na.rm = T)}, 
                  fun.max = function(x) {mean(x) + sd(x)}, 
                  geom = "errorbar", lty =1 , size =0.75, width = 0.25, position = position_dodge(width = 0.8)) +
    geom_jitter( alpha = 0.4, size = 2, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8)) +
    xlab("") +   scale_y_continuous(limits = ylim, expand = c(0,0)) + ylab(lab.y) +  scale_x_discrete(labels = barLabs) +
    scale_fill_manual(values = colours) +
    #ylim(NA,3.2)  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black", size=1),
          axis.title.x=element_blank(),
          #axis.text.x=element_blank()),
          axis.text.x=element_text(face = faces),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length.x = unit(5, "points"),
          axis.ticks.length.y = unit(5,"points") ,
          axis.title.y = element_blank(),
          axis.text = element_text(colour = "black", size=20),
          legend.position = "none",
          plot.margin = unit(c(20,0,10,0), "points")
    )
  #Full anova
  data$genotreat = paste(data$genotype,data$treatment, sep = "_")
  anovaModel = aov(get(selectColumn) ~ genotreat, data = data)
  print(HSD.test(anovaModel, alpha=0.05, "genotreat", console=F)$groups)
  
  if(graph ==T){
    exptID = readline(prompt = "Enter the experiment ID: ")
    ggsave(file = paste(selectColumn, "SAQuant", paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}
