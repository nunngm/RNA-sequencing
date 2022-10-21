## Preformed biofilm assay R script
## Made to take 96-well plate data (with plate layouts and make it into long form then graph it)

library(ggplot2)
library(tidyverse)
library(tidyr)

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
graphData = df
toi = "1" #treatment of interest
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

p = ggplot(graphData, aes(x = har_time, y = biofilm, color = treatment)) + 
  geom_jitter(data = graphData[graphData$dup<1,], width = 0.1, height = 0, alpha = 0.4) + 
  stat_summary(aes(y = biofilm, group = treat_time, colour = treatment), fun = mean, geom = 'line', size = 1) + 
  stat_summary(data = graphData[graphData$dup<1,], aes(y = biofilm, group = treat_time), fun = mean,
               fun.min = function(x) {mean(x) - sd(x)}, 
               fun.max = function(x) {mean(x) + sd(x)}, 
               geom = "errorbar", lty =1 , size =0.75, width = 0.25) +
  stat_summary(data = graphData[graphData$dup<1,], aes(y = biofilm, group = treat_time, colour = treatment), fun = mean, geom = 'point', size = 4)
p
