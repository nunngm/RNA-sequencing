##Quantification of SA levels
library(ggplot2)
library(dplyr)
library(plyr)
library(agricolae)
library(stringr)
library(emmeans)
library(investr)
library(tidyr)
set.seed(31138)

setwd()

mydata= read.table(file= "clipboard",sep= "\t",header =T) # First row is the column information where "time" is the rounded time of the measurement and well sample information is the column names for the respective data
mydata$time = as.numeric(mydata$time)
mydata[, 4:ncol(mydata)] = lapply(mydata[,4:ncol(mydata)], as.numeric)


# Blank all the columns
mydata[,4:ncol(mydata)] = lapply(mydata[,4:ncol(mydata)], function(x){ temp = x-rowMeans(mydata[,grepl("BLANK", colnames(mydata))])})
mydata = mydata[, !grepl("BLANK", colnames(mydata))] #remove blank columns
sampleGroups = as.data.frame(str_split_fixed(as.character(colnames(mydata[,4:ncol(mydata)])), pattern = "_", 3))
colnames(sampleGroups) = c("strain", "treatment", "BR")

# START FROM HERE NOW
### next steps are to convert below here is to figure out how to compare n number of nls fits and get p-values compared to wild-type control
## next steps, convert data to long format with a sample type column, group_by() sample type and make the models then.
mydata.long = gather(mydata[,3:ncol(mydata)], key = "well", value = "OD600", -time) # converts to long form
sampleGroups = as.data.frame(str_split_fixed(as.character(mydata.long[,2]), pattern = "_", 3))
colnames(sampleGroups) = c("strain", "treatment", "BR")
mydata.long = cbind(mydata.long, sampleGroups)

## testing different types of curves
testset = mydata.long[mydata.long$strain == "EV" & mydata.long$treatment =="UN",]
time = mydata.long[mydata.long$strain == "EV" & mydata.long$treatment =="UN",]$time
OD600 = mydata.long[mydata.long$strain == "EV" & mydata.long$treatment =="UN",]$OD600
ggplot() +  geom_point(aes(x=time, y=OD600),size=2) + theme_classic()

gompertz.NLL = function(a1, b1, c1, sigma) {
  size.gomp <- a1*exp(-b1*exp(-c1*time))
  -sum(dnorm(OD600, mean = size.gomp, sd = sigma, log=T))}

## get the maximum likelihood by guessing the initial values
### a1 = asymtote, b1 = x value at inflection point (~1/e up to the asymptote), c1 = rate of change, sigma= stdev
gompertz.MLE.1 <- mle2(gompertz.NLL, start = list(a1 = 0.9,b1 = 5.5, c1 = 2, n0 = 0.025, sigma = 0.05), 
                       method="L-BFGS-B", lower=0.00001)
summary(gompertz.MLE.1)

gompertz.MLE.2 <- mle2(gompertz.NLL, start = list(a1 = 0.98,b1 = 5.42, c1 = 0.085, sigma = 0.019),
                       method="BFGS")

summary(gompertz.MLE.2) #if the initial parameters were picked well summary(gompertz.MLE.1)==summary(gompertz.MLE.2)

AIC(gompertz.MLE.2)
ft.gomp <- coef(gompertz.MLE.2)
ft.gomp

prof.gomp <- profile(gompertz.MLE.2)
confint(prof.gomp)
plot(prof.gomp, 
     conf = c(99, 95, 90, 80, 50)/100, absVal=T) # profiles with confidence intervals

gompertz.NLLv2 = function(a1, b1, c1, n0, sigma) {
  size.gomp <- a1*exp(-b1*exp(-c1*time))+ n0
  -sum(dnorm(OD600, mean = size.gomp, sd = sigma, log=T))}

gompertz.MLE.v2.1 <- mle2(gompertz.NLLv2, start = list(a1 = 0.9,b1 = 5.5, c1 = 2, n0 = 0.025, sigma = 0.05), 
                       method="L-BFGS-B", lower=0.00001)
summary(gompertz.MLE.v2.1)

gompertz.MLE.v2.2 <- mle2(gompertz.NLLv2, start = list(a1 = 0.9063,b1 = 7.451, c1 = 0.09741,n0= 0.041, sigma = 0.0146),
                       method="BFGS")

summary(gompertz.MLE.v2.2) #if the initial parameters were picked well summary(gompertz.MLE.1)==summary(gompertz.MLE.2)

AIC(gompertz.MLE.v2.2)
ft.gomp.v2 <- coef(gompertz.MLE.v2.2)
ft.gomp.v2

prof.gomp.v2 <- profile(gompertz.MLE.v2.2)
confint(prof.gomp.v2)
plot(prof.gomp.v2, 
     conf = c(99, 95, 90, 80, 50)/100, absVal=T) # profiles with confidence intervals

## Trying a logistic curve
### a = shifts the curve left and right, b controls the steepness of the curve
logistic.NLL = function(a, b, sigma){
  size.logis = exp(a+(b*time))/(1+exp(a+(b*time)))
  -sum(dnorm(OD600, mean = size.logis, sd = sigma, log = T))
}
logistic.MLE.1 <- mle2(logistic.NLL, start = list(a = 3.09,b = 0.12, sigma = 0.029)
                       , method="L-BFGS-B"
                       #, lower=0.00001
                       )
summary(logistic.MLE.1)

logistic.MLE.2 <- mle2(logistic.NLL, start = list(a = -3.09,b = 0.12, sigma = 0.029),
                       method="BFGS")

summary(logistic.MLE.2) #if the initial parameters were picked well summary(logistic.MLE.1)==summary(logistic.MLE.2)

AIC(logistic.MLE.2)
ft.logis <- coef(logistic.MLE.2)
ft.logis

prof.logis <- profile(logistic.MLE.2)
confint(prof.logis)
plot(prof.logis, 
     conf = c(99, 95, 90, 80, 50)/100, absVal=T) # profiles with confidence intervals

### K = population maximum, r= rate of grwoth, n0 = population at t0
poplogis.NLL = function(K, r, n0, sigma){
  size.poplogis = K/(1+(K/n0-1)*exp(-r*time))
  -sum(dnorm(OD600, mean = size.poplogis, sd = sigma, log = T))
}

poplogis.MLE.1 <- mle2(poplogis.NLL, start = list(K = 1,r = 0.2, n0 = 0.025, sigma = 0.5)
                       , method="L-BFGS-B"
                       , lower=0.00001
)
summary(poplogis.MLE.1)

poplogis.MLE.2 <- mle2(poplogis.NLL, start = list(K = 0.89, r= 0.15, n0 = 0.027, sigma = 0.01),
                       method="BFGS")

summary(poplogis.MLE.2) #if the initial parameters were picked well summary(poplogis.MLE.1)==summary(poplogis.MLE.2)

AIC(poplogis.MLE.2)
ft.poplogis <- coef(poplogis.MLE.2)
ft.poplogis

prof.poplogis <- profile(poplogis.MLE.2)
confint(prof.poplogis)
plot(prof.poplogis, 
     conf = c(99, 95, 90, 80, 50)/100, absVal=T) # profiles with confidence intervals

# Comparing between curve types
plot(OD600 ~ time, pch = 20)
curve(ft.gomp[1]*exp(-ft.gomp[2]*exp(-ft.gomp[3]*x)), 
      from = 0, to = 48, add = T, col="red", lwd = 2)
curve(ft.gomp2[1]*exp(-ft.gomp2[2]*exp(-ft.gomp2[3]*x))+ft.gomp[4], 
      from = 0, to = 48, add = T, col="green", lwd = 2)
curve(exp(ft.logis[1]+ft.logis[2]*x)/(1+exp(ft.logis[1]+ft.logis[2]*x)), 
      from = 0, to = 48, add = T, col="blue", lwd = 2)
curve(ft.poplogis[1]/(1+(ft.poplogis[1]/ft.poplogis[3]-1)*exp(-ft.poplogis[2]*x)), 
      from = 0, to = 48, add = T, col="purple", lwd = 2)
legend("bottomright", 
       legend = c("Gompertz", "Gompertz with n0", "Logistic", "Population logistic"), 
       col = c("red", "green", "blue", "purple"), lty = 1)

AICctab(gompertz.MLE.2,gompertz.MLE.4, logistic.MLE.2, poplogis.MLE.2, nobs = 546 )
BICtab(gompertz.MLE.2,gompertz.MLE.4, logistic.MLE.2, poplogis.MLE.2, nobs = 546)

#### Examining the fit of the curves to the data the population logistic curve obviously fits the data the best.

# Let's get to the science
testset = mydata.long[mydata.long$strain == "EV",]
testset = testset[testset$treatment =="0" | testset$treatment =="0.5", ]
testset = testset[testset$treatment =="0", ]
testset$treatment = factor(testset$treatment, levels = c("0", "0.5"))

poplogis.NLL = function(K, r, n0, sigma){
  size.poplogis = K/(1+(K/n0-1)*exp(-r*time))
  -sum(dnorm(OD600, mean = size.poplogis, sd = sigma, log = T))
}


poplogis.null <- mle2(poplogis.NLL, start = list(K = 1,r = 0.2, n0 = 0.025, sigma = 0.5)
                      , method="L-BFGS-B"
                      , lower=0.00001, data = testset
)
poplogis.MLE.null <- mle2(poplogis.NLL, 
                          start=list(K = 0.57, r = 0.13, n0 = 0.024, sigma = 0.142), 
                          method="BFGS", data = testset)
summary(poplogis.MLE.null)


logLik(poplogis.MLE.null)
ft.poplogis <- coef(poplogis.MLE.null)
ft.poplogis

prof.poplogis <- profile(poplogis.MLE.null)
confint(prof.poplogis)
plot(prof.poplogis,conf = c(99, 95, 90, 80, 50)/100, absVal=T) # profiles with confidence intervals


with(testset,
     plot(OD600 ~ time, col=c("purple", "blue")[treatment], pch=16, cex=1.6))
# Females in purple, males in blue
curve(ft.poplogis[1]/(1+(ft.poplogis[1]/ft.poplogis[3]-1)*exp(-ft.poplogis[2]*x)), 
      from = 0, to = 48, add = T, col="red", lwd = 2)

  legend(x=0, y=0.8, legend=c("0", "0.5"), pch=16,col=c("purple", "blue") )


fit = nls(OD600 ~ SSfpl(time, A, B, xmid, scal), data = mydata.long[mydata.long$treatment == "gluc" & mydata.long$time%%2 ==0 ,] ) # reducing the timepoints (and outlying first bio rep really helped this model work better)
fit = nls(OD600 ~ SSfpl(time, A, B, xmid, scal), data = mydata.long[mydata.long$treatment == "0",])
#fit = nls(OD600 ~ SSlogis(time, Asym, xmid, scal), data = mydata.long[mydata.long$treatment == "0",]) # for a three parameter logistic curve but the four-parameter logistic had a better fit
#fit = nls(OD600 ~ SSgompertz(time, Asym, b2, b3), data = mydata.long[mydata.long$treatment == "0",]) # tried gompertz and it had a worse fit but probably because I wasn't setting the start conditions.
new.data = data.frame(time = seq(0, 48, length.out = 25))
interval = as_tibble(predFit(fit, newdata = new.data, interval = "confidence", level = 0.99)) %>% mutate(time = new.data$time)
 p1 <-   ggplot(mydata.long[mydata.long$treatment == "UN",]) +  geom_point(aes(x=time, y=OD600, color = strain),size=2) + xlab("Time (h)") + ylab("Optical density (OD600)") 
p1 + geom_line(data=interval, aes(x = time, y = fit ))+
  geom_ribbon(data=interval, aes(x=time, ymin=lwr, ymax=upr), alpha=0.5, inherit.aes=F, fill="blue")+
  theme_classic()

p1 + geom_line(data=new.data, aes(x = time, y = fit ))+
  geom_ribbon(data=new.data, aes(x=time, ymin=lwr, ymax=upr), alpha=0.5, inherit.aes=F, fill="blue")+
  theme_classic()

## Boot strap confidence intervals -> This is the way to go for CIs
bootFun = function(newdata){
  start = coef(fit)
  df = mydata.long[mydata.long$treatment == "gluc"& mydata.long$time%%2 ==0,]
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
new.data = data.frame(time = seq(0, 48, length.out = 25))
new.data$fit = predict(fit, newdata = new.data)
new.data$lwr <- apply(bmat, 1, quantile, 0.01, na.rm = TRUE) 
new.data$upr <- apply(bmat, 1, quantile, 0.99, na.rm = TRUE) # upr and lwr together make a confidence interval with a probability width of 0.01


# summarizing and graphing
## summarize the data in long format
dataSummary = function(data, varName, keyVariables = c("strain", "treatment", "time"), control = "0", var_interest, interval = 2){ # a simple function which just pulls out the mean optical density for each condition at the indicated intervals relative to the selected treatment

  summary_func = function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  } # collects the means
  
  data = subset(data, time%%interval==0) # only pick the times at the intervals of interest
  data_sum =ddply(data, keyVariables, .fun=summary_func,
                  varName) # summarize the data based on the key variables defined above
  data_sum = rename(data_sum, c("mean" = varName)) # rename the mean column to the original varibale's name
  data_sum = data_sum[data_sum$treatment == control|data_sum$treatment == var_interest,] # pick out only treatments of interest (e.g. the control and one other)
  
  return(data_sum)
}
data.sub = dataSummary(mydata.long, varName = "OD600", var_interest = "0.5")

ggplot(data.sub, aes(x = time, y = OD600, color = strain, shape = treatment))+geom_line() + geom_point() + geom_errorbar(aes(ymin = OD600 - sd, ymax = OD600 + sd), width = .2, position= position_dodge(0.05))

## quick bar graph at a specific time point
ggplot(data.sub[data.sub$time ==48,], aes(x=strain, y=OD600, fill = treatment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=OD600-sd, ymax=OD600+sd), width=.2,
                position=position_dodge(.9)) 

# ARETe
# this is from above you need take the "well" column and split into a strain column and a SA conc column
sampleGroups = as.data.frame(str_split_fixed(as.character(colnames(mydata[,4:ncol(mydata)])), pattern = "_", 3)) # this is from above you need take the "well" column and split into a strain column and a SA conc column

sum.bar = function(data, varName, keyVariables = c("strain", "treatment", "time"), test, control, select_time, relative = F){
  controlSummary = function(data, varName, keyVariables = c("strain", "treatment", "time"), control){
    
    summary_func = function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    
    data_sum =ddply(data, keyVariables, .fun=summary_func,
                    varName)
    data_sum = rename(data_sum, c("mean" = varName))
    data_sum = data_sum[data_sum$treatment == control,]
    
    return(data_sum)
  }
  
  ref = controlSummary(mydata.long, varName = "OD600", control = "0")
  
  data = data[data$time == select_time, ] # select the specific time point to be used
  #data = data[data$treatment == control|data$treatment == test,] # select specific tretaments
  
  if (relative ==T){
    data[data$strain == "EV", ]$OD600 = (data[data$strain == "EV", ]$OD600 - ref[ref$strain =="EV" &ref$time==0,]$OD600)/(ref[ref$strain == "EV"& ref$time == select_time,]$OD600- ref[ref$strain =="EV" &ref$time==0,]$OD600) * 100
    
    for(i in levels(as.factor(data$strain))){
      data[data$strain == i, ]$OD600 = (data[data$strain == i, ]$OD600 - ref[ref$strain ==i &ref$time==0,]$OD600)/(ref[ref$strain == i& ref$time == select_time,]$OD600- ref[ref$strain ==i &ref$time==0,]$OD600) * 100
      
    }
    xlab = "Relative growth (%)"
  }
  ### Subset by treatment
  
  ### plot the findings (under construction)
  p = ggplot(data, aes(x=strain, y=OD600, fill = treatment )) + geom_jitter( color= expCol,
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
}

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
