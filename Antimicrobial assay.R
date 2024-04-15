##Quantification of SA levels
library(ggplot2)
library(dplyr)
library(agricolae)
library(stringr)
library(emmeans)
library(investr)
library(tidyr)
set.seed(31138)

mydata= read.table(file= "clipboard",sep= "\t",header =T)
mydata$time = as.numeric(mydata$time)
mydata[, 4:ncol(mydata)] = lapply(mydata[,4:ncol(mydata)], as.numeric)

# Blank all the columns
mydata[,4:ncol(mydata)] = lapply(mydata[,4:ncol(mydata)], function(x){ temp = x-rowMeans(mydata[,grepl("BLANK", colnames(mydata))])})
mydata = mydata[, !grepl("BLANK", colnames(mydata))] #remove blank columns
sampleGroups = as.data.frame(str_split_fixed(as.character(colnames(mydata[,4:ncol(mydata)])), pattern = "_", 3))
colnames(sampleGroups) = c("strain", "AB_conc", "BR")

x = mydata$time
y = mydata$YMM12_0_1
fit <- nls(y ~ SSlogis(x, Asym, xmid, scal), trace = T, data = data.frame(x, y))

fit <- nls(mydata$YMM12_0_1 ~ SSlogis(mydata$time, Asym, xmid, scal), data = data.frame(mydata$time, mydata$YMM12_0_1))
model = lm(lm(y~I(x^3)+I(x^2)+x))
#model = lm(YMM12_0_1 ~ time, data = mydata)
mydata <- mydata %>%
  dplyr::select( -matches('fit'), -matches('lwr'), -matches('upr') ) %>%
  cbind( predict(model, interval='conf'))
ggplot(mydata, aes(x=time, y=YMM12_0_1)) +
  geom_point() 

+
  geom_ribbon( aes(ymin=lwr, ymax=upr), alpha=.3 ) +
  geom_line( aes(y=fit)) 

## The below worked



fit = nls(YMM12_500_1 ~ SSlogis(time, Asym, xmid, scal), data = mydata)
fit = nls(YMM12_500_1 ~ SSfpl(time, A, B, xmid, scal), data = mydata)
fit = nls(YMM12_0_1 ~ SSgompertz(time, Asym, b2, b3), data = mydata)

plot(mydata$YMM12_100_1~mydata$time)
lines(seq(0.5, 48, length.out = 100), predict(fit, newdata = data.frame(time = seq(0.5, 48, length.out = 100)))) # Predict just plops out the predicted y-axis points
new.data = data.frame(time = seq(0.5, 48, length.out = 100))
interval = as_tibble(predFit(fit, newdata = new.data, interval = "confidence", level = 0.95)) %>% mutate(time = new.data$time)


p1 <- ggplot(mydata) +  geom_point(aes(x=time, y=YMM12_500_1),size=2, colour="black") + xlab("Time (h)") + ylab("Optical density (OD600)") 
p1 + geom_line(data=interval, aes(x = time, y = fit ))+
  geom_ribbon(data=interval, aes(x=time, ymin=lwr, ymax=upr), alpha=1, inherit.aes=F, fill="blue")+
  theme_classic()
# START FROM HERE NOW
### next steps are to convert below here is to figure out how to compare n number of nls fits and get p-values compared to wild-type control
## next steps, convert data to long format with a sample type column, group_by() sample type and make the models then.
mydata.long = gather(mydata[,3:ncol(mydata)], key = "well", value = "OD600", -time) # converts to long form
sampleGroups = as.data.frame(str_split_fixed(as.character(mydata.long[,2]), pattern = "_", 3))
colnames(sampleGroups) = c("strain", "AB_conc", "BR")
mydata.long = cbind(mydata.long, sampleGroups)

fit = nls(OD600 ~ SSfpl(time, A, B, xmid, scal), data = mydata.long[mydata.long$AB_conc == "0" & mydata.long$time%%2 ==0 & mydata.long$BR!=1,] ) # reducing the timepoints (and outlying first bio rep really helped this model work better)
fit = nls(OD600 ~ SSfpl(time, A, B, xmid, scal), data = mydata.long[mydata.long$AB_conc == "0",])
#fit = nls(OD600 ~ SSlogis(time, Asym, xmid, scal), data = mydata.long[mydata.long$AB_conc == "0",]) # for a three parameter logistic curve but the four-parameter logistic had a better fit
#fit = nls(OD600 ~ SSgompertz(time, Asym, b2, b3), data = mydata.long[mydata.long$AB_conc == "0",]) # tried gompertz and it had a worse fit but probably because I wasn't setting the start conditions.
new.data = data.frame(time = seq(0.5, 48, length.out = 100))
interval = as_tibble(predFit(fit, newdata = new.data, interval = "confidence", level = 0.99)) %>% mutate(time = new.data$time)
p1 <- ggplot(mydata.long[mydata.long$AB_conc == "0"|mydata.long$AB_conc == "500",]) +  geom_point(aes(x=time, y=OD600),size=2, colour="black") + xlab("Time (h)") + ylab("Optical density (OD600)") 
p1 + geom_line(data=interval, aes(x = time, y = fit ))+
  geom_ribbon(data=interval, aes(x=time, ymin=lwr, ymax=upr), alpha=0.5, inherit.aes=F, fill="blue")+
  theme_classic()

p1 + geom_line(data=new.data, aes(x = time, y = fit ))+
  geom_ribbon(data=new.data, aes(x=time, ymin=lwr, ymax=upr), alpha=0.5, inherit.aes=F, fill="blue")+
  theme_classic()

## Boot strap confidence intervals -> This is the way to go for CIs
bootFun = function(newdata){
  start = coef(fit)
  df = mydata.long[mydata.long$AB_conc == "0"& mydata.long$time%%2 ==0 & mydata.long$BR!=1,]
  dfboot <- df[sample(nrow(df), size = nrow(df), replace = TRUE),]
  bootfit = try(update(fit,
                       start = start,
                       data = dfboot),
                silent = F)
  if(inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(500, bootFun(new.data)) # run the boot strapping
### Make the line to graph on your plot
new.data = data.frame(time = seq(0.5, 48, length.out = 100))
new.data$fit = new$fit <- predict(fit, newdata = new.data)
new.data$lwr <- apply(bmat, 1, quantile, 0.005, na.rm = TRUE) 
new.data$upr <- apply(bmat, 1, quantile, 0.995, na.rm = TRUE) # upr and lwr together make a confidence interval with a probability width of 0.01


# ARETE
# this is from above you need take the "well" column and split into a strain column and a SA conc column
sampleGroups = as.data.frame(str_split_fixed(as.character(colnames(mydata[,4:ncol(mydata)])), pattern = "_", 3)) # this is from above you need take the "well" column and split into a strain column and a SA conc column


toGraph = t(mydata[, 4:ncol(mydata)])
colnames(toGraph) = mydata$time

SAQuantGraph = function(data, genotypeCol = c("#378717","#FFFF00", "#6DFDFD"), ylim = c(4,8), expCol = NA, graph = F, width = 5, height = 4, exptID = "temp", box = F){
}
SAQuantGraph = function(data,genotypes = c("Col-0", "ald1-T2", "fmo1-1"), colours = c("#378717","#FFFF00", "#6DFDFD"), selectTimes = c("UN", "12", "18", "24"), 
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
SAQuantGraph(temp, genotypes = c("Col-0", "pen3-4", "pdr12-3", "p3p12"), colours = c("#FFFFFF", "#00BBFF", "#FF3853", "#5D2169"), selectTimes = c("UN", "24"), barLabs = c("Untreated", "24 hpi"), SA = "inter", ylim = c(0,2500), graph = T, height = 6, width = 8)

#Double mutant colours
colours = c("#FFFFFF", "#00BBFF", "#FF3853", "#5D2169")
