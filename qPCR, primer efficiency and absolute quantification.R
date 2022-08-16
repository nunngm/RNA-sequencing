## Absolute quantification
### Based on Rutledge and Stewart, 2008
# Optical calibration factor (FU/ng dsDNA) for our kit (Biorad Universal SYBR green) = 957 FU/ng dsDNA
# F0 = calculated initial fluorescence, M0 = calculated initial mass of target double-stranded cDNA, N0 number of initial transcripts, OCF = optical calibration factor, As = amplicon size
# The N0 is scaled by amplicon size irrespective of GC content
amplicon_sizes =c(ALD1 = 105,FMO1 = 124, SEC5A = 93, CUL4 = 106,ACTIN1 = 169, UGT76B1 = 121,SARD1 = 139, CBP60G = 101, FLS2 = 112, ICS1 = 154, SOC1 = 136, SVP = 122, UGT74F1 = 116, PR1 = 89, RLP28 = 119)

library(qpcR)
library(tidyr)
library(ggplot2)
library(dplyr)
library(svglite)
library(agricolae)
library(car) 

mydata= read.table(file= "clipboard",sep= "\t",header =T)
samples = as.character(mydata[1,2:ncol(mydata)])
mydata = mydata[2:nrow(mydata),]
mydata = as.data.frame(lapply(mydata, as.numeric))
groups = as.data.frame(t(as.data.frame(strsplit(samples, split = "_"))))
colnames(groups) = c("gene", "type","rep")
sample_group = paste(groups$gene, groups$type, sep = "_")

#laptop directory
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-qRT-PCR\\Primer Efficiency test")

#Desktop directory
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-qRT-PCR\\Primer Efficiency test")

curves = lapply(2:ncol(mydata), function(x){
     pcrfit(mydata, cyc = 1, fluo = x)})
## Determine the average threshold cycle
E = data.frame(t(sapply(1:(ncol(mydata)-1), simplify = "array", function(x){efficiency(curves[[x]], type = "cpD2")})))
unlist(E$eff)/2 #
Ft = mean(unlist(E$fluo[E$fluo>100])) ##Might error out do to wells that don't amplify
selSamples = mydata[40,2:ncol(mydata)]>Ft
groups = groups[selSamples,]
mydata = mydata[,c(T, selSamples)]
curves = lapply(2:ncol(mydata), function(x){pcrfit(mydata, cyc = 1, fluo = x)})

##Using the set threshold cycle find the cycle in which each well passes the threshold
E = data.frame(t(sapply(1:(ncol(mydata)-1), simplify = "array", function(x){efficiency(curves[[x]], type = "cpD2", thresh = Ft)})))
Ct = unlist(E$cpT)
Emax = unlist(E$eff)/2 ##Generates the amplification frequencies of each well

## Check the efficicies and build the df for doing differential expression analysis
df = cbind(groups, Emax, Ct)
aggregate(Emax ~ gene, df, mean) #check the means for each of the genes, should be 1+-0.1
aggregate(Ct ~ gene, df, mean)
df = df[ order(df$gene, df$type, df$rep), ]#The data frame MUST be ordered correctly before starting this analysis all samples need to be in the correct location, check they are match before running
df = cbind(df,sample = paste(df[df$gene == df$gene[], ]$type, df[df$gene == df$gene[], ]$rep, sep = "_"))

# convert to wide format
mat = matrix(NA, ncol = length(levels(as.factor(df$gene))), nrow = length(paste(df[df$gene == refGeneName, ]$type, df[df$gene == refGeneName, ]$rep, sep = "_")))
colnames(mat) = levels(as.factor(df$gene))
rownames(mat) = paste(df[df$gene == refGeneName, ]$type, df[df$gene == refGeneName, ]$rep, sep = "_")
for (i in 1:nrow(df)){
     mat[rownames(mat) == df$sample[i] ,colnames(mat)==df$gene[i]] = df$Ct[i]
}
df.wide = as.data.frame(mat)

refGeneName = "SEC5A"
targetGeneName = "RLP28"

#Differential expression analysis 
refGeneExp = setNames(df.wide[[refGeneName]], rownames(df.wide)) #Change this to enter your reference gene expression (if from another plate that worked well you will be able to compare to the whole plate)
                        
df.wide = cbind(df.wide, refGeneExp)
ampSize = c(amplicon_sizes[colnames(df.wide[, 1:ncol(df.wide)-1])], amplicon_sizes[refGeneName]) #binds the ampSize for each of the products and the reference gene (need to change after the comma if using a different reference gene)


df.wide = sapply(df.wide[, 1:ncol(df.wide)], as.numeric, simplify = "array")
df.wide = as.data.frame(t(log2(2^(t(df.wide[, 1:ncol(df.wide)])*-1) / ampSize * amplicon_sizes[refGeneName])* -1)) # scale expression of expre based on the amplicon size of the reference gene
rownames(df.wide) = paste(df[df$gene == refGeneName, ]$type, df[df$gene == refGeneName, ]$rep, sep = "_")
# df.wide = as.data.frame(cbind(df.wide, type = df[df$gene == df$gene[1], ]$type))

expression = cbind(as.data.frame(lapply(df.wide[,1:ncol(df.wide)-1], function(x){2^-(as.numeric(x)-as.numeric(df.wide$refGeneExp))}), row.names = rownames(df.wide)), t(as.data.frame(strsplit(rownames(df.wide), split = "_"), row.names = c("type", "rep"))), stringsAsFactors = T) # Takes the wide-format data frame and calculates the expression for the genes relative to the reference gene.
aggregate(. ~ type, expression[, 1:ncol(expression)-1], mean) # quick check the means to see expression

# Clean up the data frame
expression$type = as.factor(as.integer(expression$type)+3)
#expression = expression[7:nrow(expression),] # removing 1-2 wpg samples

## Plotting
is_outlier <- function(x) {
  return(x < quantile(x, 0.25, na.rm = T) - 1.5 * IQR(x, na.rm=T) | x > quantile(x, 0.70, na.rm = T) + 1.5 * IQR(x,na.rm = T))
} ## Function for determining outliers
     
p = expression %>%
     group_by(type) %>%
     mutate(inlier = ifelse(is_outlier(!!as.name(targetGeneName)), as.numeric(NA), !!as.name(targetGeneName)), outlier = ifelse(is_outlier(!!as.name(targetGeneName)), !!as.name(targetGeneName), as.numeric(NA)) ) %>%
     ggplot(., aes(x=type, y=inlier, colour = rep)) +
     stat_summary(fun = mean, geom = "bar", fill = c( #"#444444",
          "#666666", "#9A9A9A", "#CDCDCD", "#FFFFFF"), colour = "#000000", size = 0.75) +
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

## Stats
# Label inliers
analysis = expression %>%
     group_by(type) %>%
     mutate(inlier = ifelse(is_outlier(!!as.name(targetGeneName)), as.numeric(NA), !!as.name(targetGeneName)))
analysis = analysis[is.na(analysis$inlier)==F,] #remove outliers

# label each row with the number of reps
temp = summary(analysis)
temp = as.data.frame(t(data.frame(strsplit(temp[,"type"], split = ":", fixed = T))))
temp = setNames(as.integer(temp[, 2]),temp[,1])
analysis = cbind(analysis,reps =temp[analysis$type])
rm(temp)

#Automated
AutomatedStats = function(md, noTukey = F){ 
  results = matrix(NA,ncol =4,nrow =2)
  colnames(results) = c("ANOVA","Shapiro-Wilk","Levene","Bartlett")
  rownames(results) = c("P-value","Result")
  model = paste0(targetGeneName, " ~ type")
  anovaModel = aov(as.formula(model), data=md) #the numeric variable must be first in this expression 
  r = summary(anovaModel)
  results[1,1] = r[[1]]$`Pr(>F)`[1]
  md.lm = lm(as.formula(model), data = md)
  md.res = resid(md.lm)
  results[1,2]=shapiro.test(md.res)[2]$p.value
  library(car) 
  results[1,3]=leveneTest(as.formula(model), data=md)[3][1,1]
  results[1,4] = bartlett.test(as.formula(model), data=md)[3]$p.value
  par(mfrow = c(2,2))
  #plot(md.lm)
  par(mfrow = c(1,1))
  library(agricolae)
  results[2,] = rep("FAIL",each = ncol(results))
  if (as.numeric(results[1,1])<=0.05){
    results[2,1]= "PASS"
  }
  for (i in 2:ncol(results)){
    if (as.numeric(results[1,i])>0.05){
      results[2,i]="PASS"
    }
  }
  results[1,]= signif(as.double(results[1,]), digits=4)
  if (noTukey ==T){
    return(results)
  }
  y=NULL #this was designed to be an output variable if you need to output both stats table and tukey grouping
  y$groups= HSD.test(anovaModel, alpha=0.05, "type", console=T)$groups#gives letter codes but no p values
  y$stats=results
  return(y)
}
stats = AutomatedStats(analysis)
analysis[[targetGeneName]] = log2(analysis[[targetGeneName]])


genesOI = levels(as.factor(df$gene)) #Make sure to keep the reference gene in so it is checked every time.
lapply(genesOI, simplify = "array", function(x){
     df[df$gene == x, ]$Ct - rE$Ct
})

paste0 

## Log graph
temp = expression
temp[[targetGeneName]] = log2(temp[[targetGeneName]]) # log2 transform
m <- floor(min(temp[[targetGeneName]])) # find min
temp[[targetGeneName]] <- temp[[targetGeneName]] - m

##plotting
p = temp %>%
     group_by(type) %>%
     mutate(inlier = ifelse(is_outlier(!!as.name(targetGeneName)), as.numeric(NA), !!as.name(targetGeneName)), outlier = ifelse(is_outlier(!!as.name(targetGeneName)), !!as.name(targetGeneName), as.numeric(NA)) ) %>%
     ggplot(., aes(x=type, y=inlier, colour = rep)) +
     stat_summary(fun = mean, geom = "bar", fill = c( "#444444",
          "#666666", "#9A9A9A", "#CDCDCD", "#FFFFFF"), colour = "#000000", size = 0.75) +
     stat_summary(fun = mean,
          fun.min = function(x) {mean(x) - sd(x)}, 
          fun.max = function(x) {mean(x) + sd(x)}, 
          geom = "errorbar", lty =1 , size =0.75, width = 0.25, colour = "#000000") +
     #geom_boxplot(fill = rep(c("#FFFFFF"), 5)) +
     geom_jitter(width = 0.25, color= "#000000", size = 2, alpha = 0.4) + scale_y_continuous(expand = expansion(c(0, 0.1)), breaks=seq(floor(min(temp[[targetGeneName]])), ceiling(max(temp[[targetGeneName]])), length=5),labels=as.character(round(seq(m, ceiling(max(temp[[targetGeneName]]+m)), length=5),2))) + 
     geom_point(aes(x = type, y = outlier), size =2, alpha = 1, shape = 8, colour = "#000000")   +
theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) + 
    xlab("Weeks post-germination (wpg)") + ylab(paste0(targetGeneName,"/SEC5A (log2)")) + theme(panel.grid.major = element_blank(),  
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
ggsave(file = paste0(targetGeneName,"_WKEX-22-1-log.svg"), plot = p, width = 5, height = 4)
#__________absolute quatitation_________________
##ss

E = efficiency(curves, type = "cpD2",thresh = Ft)
Ct = E$cpT
Emax = E$eff/2
Ft = 187
F0 = Ft/(Emax+1)^Ct
OCF = 957
M0 = F0/OCF
N0 = M0 * 9.1 * 10^11 / As # where 9.1 x 10^11 is the number of base pairs per nanogram of dsDNA.

OCF = 356 # OCF of a similar 
As = c(rep(105, times = 24), rep(124, times = 24), rep(106, times = 24), rep(93, times = 24))

curves = lapply(2:ncol(mydata), function(x){pcrfit(mydata, cyc = 1, fluo = x)})


E = data.frame(t(sapply(1:(ncol(mydata)-1), simplify = "array", function(x){efficiency(curves[[x]], type = "cpD2")})))
Ft = mean(unlist(E$fluo))
E = data.frame(t(sapply(1:(ncol(mydata)-1), simplify = "array", function(x){efficiency(curves[[x]], type = "cpD2", thresh = Ft)})))
Ct = unlist(E$cpT)
Emax = unlist(E$eff)/2
F0 = Ft/(Emax+1)^Ct
M0 = F0/OCF
N0 = M0 * 9.1 * 10^11 / As
temp = cbind(sample_group,N0)
write.xlsx2(temp, file = "qPCR absolute abundances.xlsx")


N0 = as.data.frame(cbind(sample_group = paste(groups$gene, groups$type, sep = "_"), rep = groups$rep, N0))
N0 = reshape(N0, idvar = "rep", timevar = "sample_group",direction = "wide")
N0 = N0[,2:ncol(N0)]

colMeans(N0, na.rm = T)
summary(N0)
