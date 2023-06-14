## Preformed biofilm assay R script
## Made to take 96-well plate data (with plate layouts and make it into long form then graph it)

library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
setwd("C:\\Users\\garre\\OneDrive\\Desktop\\stuff")

#Read in the blanked plate data
biofilmData = unlist(read.table(file= "clipboard",sep= "\t",header =F)) #unlist to linearize the data

# Read in plate layout
# The information must be spatially correct (all thre wells need to have a corresponding well on the plate data)
plateLayout = unlist(read.table(file= "clipboard",sep= "\t",header =F)) #unlist to linearize the data
names(plateLayout) == names(biofilmData) #Check that you have the same number of rows/columns and wells in both data sets, does not check whether the layout is correct, that is up to you.

# Read growth data (recorded prior to sacrificing wells for growth/death)
growthData = unlist(read.table(file= "clipboard",sep= "\t",header =F))


linData = as.data.frame(cbind(names(plateLayout), plateLayout, biofilmData, growthData))
colnames(linData) = c("id", "group", "biofilm", "growth")
linData = linData[linData$group != "",] #wells which don't have group information are removed because it is assumed they are blank/mistake was made

# Make a data frame with group information split apart
df = data.frame(cbind(t(data.frame(strsplit(linData[,2], split = "_"))), linData[,3], linData[,4]))
colnames(df) = c("har_time", "treat_time", "treatment","har_type" ,"biofilm", "growth") # har_time = hours after start of experiment that the well was harvested, treat_time = the time after the start that the well was treated, treatment = concentration of SA the well was treated with (UN = untreated), har_type = what the well was used for (H = well was tested for planktonic growth; BS = after dumping&rinsing biofilm was scraped to test for growth of biofilm), OD = optical density observed in that well.
rownames(df) = linData[,1]
### BIOFILM STAINING ANALYSIS
df = df[df$har_type == "BS",] #remove any wells that were sacrificed or used for other purposes than biofilm staining

toi = "0.1" #treatment of interest
graphData = df
graphData = graphData[df$treatment == "UN"| df$treatment == toi, ]
graphData$treat_time = paste(graphData$treat_time, graphData$treatment, sep = "_")
graphData = cbind(graphData, 0)
colnames(graphData)[7] = "dup"

# Duplicate untreated data to attach lines to last untreated datapoint , could remove
dupData = graphData[graphData$treatment == "UN",]
dupData$treat_time = paste(dupData$har_time, toi, sep = "_")
dupData$treatment = toi
dupData$dup = 1

graphData = rbind(graphData, dupData)
graphData$biofilm = as.numeric(graphData$biofilm)
graphData$growth = as.numeric(graphData$growth)
graphData$treat_time = factor(graphData$treat_time, levels = c("0_UN", paste0("0_", toi), paste0("24_", toi),paste0("48_", toi), paste0("72_", toi), paste0("96_", toi) ))

graphColours = c("#E84855", "#000000", "#E84855", "#FF9B71", "#FFFD82", "#000000")
names(graphColours) = levels(graphData$treatment)

##My favourite version of the graph #line graph where each treatment time branches off from untreated
p = ggplot(graphData, aes(x = har_time, y = biofilm, colour = treatment)) + 
  geom_jitter(data = graphData[graphData$dup<1,], width = 0.05, height = 0, alpha = 0.4, size = 2) + 
  stat_summary(aes(y = biofilm, group = treat_time, colour = treatment), fun = mean, geom = 'line', size = 1) + 
  stat_summary(data = graphData[graphData$dup<1,], aes(y = biofilm, group = treat_time), fun = mean,
               fun.min = function(x) {mean(x) - sd(x)}, 
               fun.max = function(x) {mean(x) + sd(x)}, 
               geom = "errorbar", lty =1 , size =0.75, width = 0.05) +
  stat_summary(data = graphData[graphData$dup<1,], aes(y = biofilm, group = treat_time, colour = treatment), fun = mean, geom = 'point', size = 3)  +
  xlab("Time (hpi)") + ylab(bquote('Biofilm formation (OD'[570]~')')) + 
  ylim(NA,3.2)  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size=1),
        axis.title.x=element_text(size=30),
        #axis.text.x=element_blank()),
        axis.ticks=element_line(colour = "black", size =1),
        axis.ticks.length = unit(10,"points") ,
        axis.title.y = element_text(size=30),
        axis.text = element_text(colour = "black", size=30)
  ) + scale_colour_manual(values = graphColours)
p
# Bar graph showing same data as line graph but as a bar
graphData = graphData
p = ggplot(graphData, aes(x = har_time, y = biofilm, fill = treat_time)) + 
  stat_summary(aes(y = biofilm, group = treat_time, fill = treat_time), colour = "#000000", fun = mean, geom = 'bar', width = 0.75, size = 1, position = position_dodge(width = 0.75)) + 
  stat_summary( aes(y = biofilm, group = treat_time), fun = mean,
                fun.min = function(x) {mean(x) - sd(x)}, 
                fun.max = function(x) {mean(x) + sd(x)}, 
                geom = "errorbar", lty =1 , size =0.75, width = 0.05, position = position_dodge(width = 0.75)) +
  geom_jitter( alpha = 0.4, size = 2, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  xlab("Time (hpi)") + scale_y_continuous(limits = c(0,3.2), expand = c(0,0)) + ylab(bquote('Biofilm formation (OD'[570]~')')) + 
  scale_x_discrete(labels = barLabs) +
  #ylim(NA,3.2)  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size=1),
        axis.title.x=element_text(size=30),
        #axis.text.x=element_blank()),
        axis.ticks=element_line(colour = "black", size =1),
        axis.ticks.length = unit(10,"points") ,
        axis.title.y = element_text(size=30),
        axis.text = element_text(colour = "black", size=30)
  )
p
#p
svg(filename = paste0("PBA-22-1-", toi, "mM SA-bifilms-oneCol.svg"), width = 10, height = 8)
p
dev.off()

graphData = graphData[graphData$treat_time != paste0("96_",toi) & graphData$har_time !=0,] #remove 96_untreated data
barLabs = c("Untreated", "0","24", "48", "72")
p = ggplot(graphData, aes(x = treat_time, y = biofilm, fill = har_time)) + 
  stat_summary(aes(y = biofilm, group = har_time, fill = har_time), colour = "#000000", fun = mean, geom = 'bar', width = 0.75, size = 1, position = position_dodge2(width = 0.75, preserve = "single")) + 
  stat_summary( aes(y = biofilm, group = har_time), fun = mean,
               fun.min = function(x) {mean(x) - sd(x)}, 
               fun.max = function(x) {mean(x) + sd(x)}, 
               geom = "errorbar", lty =1 , size =0.75, width = 0.05, position = position_dodge(width = 0.75)) +
  geom_jitter( alpha = 0.4, size = 2, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  xlab("Time of treatment (hpi)") + scale_y_continuous(limits = c(0,3.2), expand = c(0,0)) + ylab(bquote('Biofilm formation (OD'[570]~')')) + 
  scale_x_discrete(labels = barLabs) +
  #ylim(NA,3.2)  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size=1),
        axis.title.x=element_text(size=30),
        #axis.text.x=element_blank()),
        axis.ticks=element_line(colour = "black", size =1),
        axis.ticks.length = unit(10,"points") ,
        axis.title.y = element_text(size=30),
        axis.text = element_text(colour = "black", size=30)
  )
p

ToTGraph = function(data, treatmentOfInterest = 0, graphColours = c("#FFFC54","#FF9B71", "#E84855","#962B2B") , width = 7, height = 5, biofilm = T, barLabs = c("Untreated", "0","24", "48", "72"), graph = F){
  
  if (biofilm == T){
    selectColumn = "biofilm"
    lab.y = bquote('Biofilm formation (OD'[570]~')')
  }else{
    selectColumn = "growth"
    lab.y = bquote('Bacterial growth (OD'[600]~')')
  }
  #take full data and pare down to just the required treatment and untreated
  df.select = df[df$har_type == "BS",] #remove any wells that were sacrificed or used for other purposes than biofilm staining
  df.select = df.select[df.select$treatment == "UN"| df.select$treatment == treatmentOfInterest, ]
  df.select$treat_time = paste(df.select$treat_time, df.select$treatment, sep = "_")
  df.select = df.select[df.select$treat_time != paste0("72_", treatmentOfInterest) | df.select$har_time !=96,]
  df.select = df.select[df.select$har_time != 0,]
  df.select$biofilm = as.numeric(df.select$biofilm)
  df.select$growth = as.numeric(df.select$growth)
  df.select$treat_time = factor(df.select$treat_time, levels = c("0_UN", paste0("0_", treatmentOfInterest), paste0("24_", treatmentOfInterest),paste0("48_", treatmentOfInterest), paste0("72_", treatmentOfInterest), paste0("96_", treatmentOfInterest) ))
  
  #Summarise data
  df.summary = df.select %>% group_by(treat_time, har_time) %>% summarise(mean = mean(get(selectColumn)), sd = sd(get(selectColumn)))
  
  p = ggplot(df.select, aes(x = treat_time, y = get(selectColumn), fill = har_time)) + 
    geom_bar(mapping = aes(x = treat_time, y = mean, fill = har_time), data = df.summary, stat = "identity", position = position_dodge2(width = 0.75,preserve = "total")) +
    geom_errorbar(mapping = aes(x = treat_time,y = mean, ymin = mean-sd, ymax = mean+sd), data = df.summary, size = 0.5, position = position_dodge2(width = 0.75, padding = 0.65, preserve = "total")) +
    geom_point(alpha = 0.4, size = 2, stat = "identity", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1)) +
    xlab("Time of treatment (hpi)") + scale_y_continuous(limits = c(0,3.2), expand = c(0,0)) + ylab(lab.y) + 
    scale_x_discrete(labels = barLabs) + scale_fill_manual(values = graphColours)+
    #ylim(NA,3.2)  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black", size=1),
          #axis.title.x=element_text(size=15),
          axis.title.x = element_blank(),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length = unit(5,"points") ,
          axis.title.y = element_text(size=15),
          axis.text = element_text(colour = "black", size=15)
    )
  if(graph ==T){
    exptID = readline(prompt = "Enter the experiment ID: ")
    ggsave(file = paste(selectColumn, treatmentOfInterest, paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}
## OG format of this data
toi = "0.1" #treatment time of interest
graphData = df
graphData = graphData[df$treat_time == "0" & df$treatment != "UN", ]
graphData$biofilm = as.numeric(graphData$biofilm)
graphData$growth = as.numeric(graphData$growth)
graphData$treatment = factor(graphData$treatment, levels = c("0", "0.005", "0.05", "0.1", "1", "5", "10"))

barLabs = c("0", "24","48", "72", "96")
p = ggplot(graphData, aes(x = har_time, y = biofilm, fill = treatment)) + 
  stat_summary(aes(y = biofilm, group = treatment, fill = treatment), colour = "#000000", fun = mean, geom = 'bar', width = 0.75, size = 1, position = position_dodge(width = 0.75)) + 
  stat_summary( aes(y = biofilm, group = treatment), fun = mean,
                fun.min = function(x) {mean(x) - sd(x)}, 
                fun.max = function(x) {mean(x) + sd(x)}, 
                geom = "errorbar", lty =1 , size =0.75, width = 0.05, position = position_dodge(width = 0.75)) +
  geom_jitter( alpha = 0.4, size = 2, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  xlab("Time (hpi)") + scale_y_continuous(limits = c(0,3.2), expand = c(0,0)) + ylab(bquote('Biofilm formation (OD'[570]~')')) + 
  scale_x_discrete(labels = barLabs) +
  #ylim(NA,3.2)  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size=1),
        axis.title.x=element_text(size=30),
        #axis.text.x=element_blank()),
        axis.ticks=element_line(colour = "black", size =1),
        axis.ticks.length = unit(10,"points") ,
        axis.title.y = element_text(size=30),
        axis.text = element_text(colour = "black", size=30)
  )
p

PBABarGraph = function(data, selectTreatments = c("0","0.05","0.1","1"), selectTreatTimes = c("0", "24", "48", "72"), 
                       barLabs = c("24", "48", "72", "96", "48", "72", "96", "72", "96", "96"),
                       graphColours = c("#FFFC54","#FF9B71", "#E84855","#962B2B"), width = 12, height = 5, biofilm = T, 
                       graph = F, ylim = c(0,NA)){
  if (biofilm == T){
    selectColumn = "biofilm"
    lab.y = bquote('Biofilm formation (OD'[570]~')')
  }else{
    selectColumn = "growth"
    lab.y = bquote('Bacterial growth (OD'[600]~')')
  }
  #take full data and pare down to just the required treatment and untreated
  df.select = data[data$har_type == "BS",] #remove any wells that were sacrificed or used for other purposes than biofilm staining
  df.select = df.select[as.logical(rowSums(as.data.frame(lapply(selectTreatments, function(x){df.select$treatment == x})))),]
  df.select = df.select[as.logical(rowSums(as.data.frame(lapply(selectTreatTimes, function(x){df.select$treat_time == x})))),]
  df.select$har_type = paste(df.select$treat_time, df.select$har_time, sep = "_")
  df.select$har_type = as.factor(df.select$har_type)
  df.select$biofilm = as.numeric(df.select$biofilm)
  df.select$growth = as.numeric(df.select$growth)
  
  p = ggplot(df.select, aes(x = har_type, y = get(selectColumn), fill = treatment, group = treatment)) + 
    stat_summary(aes(group = treatment), colour = "#000000", fun = mean, geom = 'bar', width = 0.75, size = 1, position = position_dodge(width = 0.75)) + 
    stat_summary( aes(y = biofilm, group = treatment), fun = mean,
                  fun.min = function(x) {mean(x) - sd(x)}, 
                  fun.max = function(x) {mean(x) + sd(x)}, 
                  geom = "errorbar", lty =1 , size =0.75, width = 0.25, position = position_dodge(width = 0.75)) +
    geom_jitter( alpha = 0.4, size = 2, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
    xlab("Time (hpi)") +   scale_y_continuous(limits = ylim, expand = c(0,0)) + ylab(lab.y) +  scale_x_discrete(labels = barLabs) +
    scale_fill_manual(values = graphColours) +
    #ylim(NA,3.2)  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black", size=1),
          axis.title.x=element_blank(),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length.x = unit(0, "points"),
          axis.ticks.length.y = unit(5,"points") ,
          axis.title.y = element_text(size=15),
          axis.text = element_text(colour = "black", size=15),
          axis.text.x = element_text(vjust = -0.25)
    )
  if(graph ==T){
    exptID = readline(prompt = "Enter the experiment ID: ")
    ggsave(file = paste(selectColumn, "PBABarGraph", paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}
