##Quantification of SA levels
library(ggplot2)
library(dplyr)
library(plyr)
library(agricolae)
library(stringr)
library(emmeans)
library(investr)
library(tidyr)
detach(package:plyr)
set.seed(31138)

setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR-NPR&PEN3\\exp-AMA\\graphs")

mydata= read.table(file= "clipboard",sep= "\t",header =T) # First row is the column information where "time" is the rounded time of the measurement and well sample information is the column names for the respective data
mydata[, 1:ncol(mydata)] = lapply(mydata[,1:ncol(mydata)], as.numeric)

# Blank all the columns
mydata[,2:ncol(mydata)] = lapply(mydata[,2:ncol(mydata)], function(x){ temp = x-rowMeans(mydata[,grepl("BLANK", colnames(mydata))])})
mydata = mydata[, !grepl("BLANK", colnames(mydata))] #remove blank columns
sampleGroups = as.data.frame(str_split_fixed(as.character(colnames(mydata[,4:ncol(mydata)])), pattern = "_", 3))
colnames(sampleGroups) = c("strain", "treatment", "BR")

# START FROM HERE NOW
### next steps are to convert below here is to figure out how to compare n number of nls fits and get p-values compared to wild-type control
## next steps, convert data to long format with a sample type column, group_by() sample type and make the models then.
mydata.long = gather(mydata[,1:ncol(mydata)], key = "well", value = "OD600", -time) # converts to long form
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

## plotting the curves
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

## Looking at the log likelihood of each model
logLik(gompertz.MLE.2)
logLik(gompertz.MLE.4)
logLik(logistic.MLE.2)
logLik(poplogis.MLE.2)

### in interpreting log likelihood results, the actual number doesn't really matter but can range from positive to negative infinity and describes how well your model fits the data it was modelled to.
### Examining the relative log likelihoods of the different models, the model that best fits the data will be the model with the highest log likelihood.

AICctab(gompertz.MLE.2,gompertz.MLE.4, logistic.MLE.2, poplogis.MLE.2, nobs = 546 )
BICtab(gompertz.MLE.2,gompertz.MLE.4, logistic.MLE.2, poplogis.MLE.2, nobs = nrow(testset))

#### Examining the fit of the curves to the data the population logistic curve obviously fits the data the best.

# Let's get to the science

## Begin with building a model of the null hypothesis, that growth is unaffected by SA concentration
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

## let's beging by using mle to test the hypothesis that only the final carrying capacity (K) is affected by SA treatment


poplogis.MLE.K = mle2(OD600 ~ dnorm(K/(1+(K/n0-1)*exp(-r*time)), sd = sigma),
                      start = list(K = ft.poplogis[1], r = ft.poplogis[2], n0 = ft.poplogis[3], sigma = ft.poplogis[4]),
                      parameters = list(K ~ treatment), method = "BFGS", data = testset
                      )
summary(poplogis.MLE.K)
ft.poplogis.K = coef(poplogis.MLE.K)

#How does the model compare to the null model
logLik(poplogis.MLE.K)
logLik(poplogis.MLE.null)

anova(poplogis.MLE.K, poplogis.MLE.null)
AICtab(poplogis.MLE.K,poplogis.MLE.null, delta=T, base=T, sort=T, weights=T)
BICtab(gompertz.MLE.a,gompertz.MLE.null, delta=T, base=T, sort=T, weights=T, nobs=nrow(hyena.growth))

with(testset,
     plot(OD600 ~ time, col=c("purple", "blue")[treatment], pch=16, cex=1))
curve(ft.poplogis[1]/(1+(ft.poplogis[1]/ft.poplogis[3]-1)*exp(-ft.poplogis[2]*x)), 
      from = 0, to = 48, add = T, col="red", lwd = 2)
curve(ft.poplogis.K[1]/(1+(ft.poplogis.K[1]/ft.poplogis.K[4]-1)*exp(-ft.poplogis.K[3]*x)), 
      from = 0, to = 48, add = T, col="green", lwd = 2)
curve((ft.poplogis.K[1]+ ft.poplogis.K[2])/(1+((ft.poplogis.K[1]+ft.poplogis.K[2])/ft.poplogis.K[4]-1)*exp(-ft.poplogis.K[3]*x)), 
      from = 0, to = 48, add = T, col="orange", lwd = 2)

### There is a difference based on carry capacity but does SA also effect the rate of growth of yeast cells


poplogis.MLE.full <- mle2(OD600 ~ dnorm(K/(1+(K/n0-1)*exp(-r*time)), sd = sigma),
                                           start = list(K = 0.8, r = 0.14, n0 = 0.025, sigma = 0.012),
                                           parameters = list(K ~ treatment, r ~ treatment), method = "BFGS", data = testset
                                    )
  
summary(poplogis.MLE.full) ## it looks like growth rate, r is also affected by treatment

prof.full <- profile(poplogis.MLE.full)
plot(prof.full)
ft.poplogis.full = coef(poplogis.MLE.full)

  mle2(BodyLength ~ dnorm(a*exp(-b*exp(-c*Age)), sd=sigma), 
                          start=list(a=97,b=0.89, c=0.11, sigma=5.57), 
                          parameters=list(a ~ Sex, b ~ Sex, c~Sex), method="BFGS", data=hyena.growth)


# taking the data in long format and graphing it

YAMA
# Below this line is not working garbage
# ---------------------------------------------------------------------------

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

## lets graph it
AMAgrowthCurve = function(df, # a data frame in long format
                          treatment, #treatments you are interested
                          strain, # strains you are interested in
                          colours = NA, # the specified colours, can be a list
                          interval = 2, 
                          width = 5, ylim = c(0,NA), height = 5, graph = F){ # graphing information
  if (length(colours) != length(strain)){
    warning("You did not include the correct number of colours, resorting to default colour scheme")
    colours = 1:length(strain)
  }
  
  df = df[df$treatment %in% treatment,]
  df$treatment = factor(df$treatment, levels = treatment)
  df = df[df$strain %in% strain,]
  df$strain = factor(df$strain, levels = strain)
  df$strainTreat = paste0(df$strain, df$treatment)
  df = df[df$time %% 2 ==0,]
  df.summary = df %>% group_by(strain, treatment, time,strainTreat) %>% summarise(growth = mean(OD600), stdev = sd(OD600))
  p = ggplot(data = df.summary, aes(x = time, y = growth, group = strainTreat)) + geom_line(aes(linetype = treatment, colour = strain),size = 1)+
    geom_errorbar(aes(ymin=growth-stdev, ymax=growth+stdev, colour = strain), width=.1, size = 1) +
    geom_point(aes(colour = strain, shape =treatment),size = 4, fill = "#FFFFFF") +
    xlab("Time (h)") + ylab(bquote('Optical density (600 nm)')) +
    scale_y_continuous(limits = ylim, expand = c(0,0)) +
    scale_x_continuous( breaks = seq(0,max(df$time, na.rm = T), by = 12)) +
    # scale_fill_manual(values = colours) +
    scale_color_manual(values = colours) +
    scale_shape_manual(values = c(16,21)) + 
    # scale_linetype_manual(values = c(1,2,2)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size=1),
          axis.title.x=element_blank(),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length.x = unit(5, "points"),
          axis.ticks.length.y = unit(5,"points") ,
          axis.title.y = element_text(size = 20),
          axis.text = element_text(colour = "black", size=15),
          axis.text.x = element_text(vjust = -0.25),
          legend.position = "none",
          plot.margin = unit(c(20,0,10,0), "points")
    )
  if(graph ==T){
    exptID = readline(prompt = "Enter the experiment ID: ")
    ggsave(file = paste("AMAgrowthCurve",paste(treatment, collapse = "and"), paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}

AMAgrowthCurve(pdr.24.3, treatment = c("0", "gluc"), ylim = c(0,1), 
               strain = c("EV", "PEN3", "PDR12"),width = 6, height = 4, colours = c("#000000", "#FF3853", "#00BBFF"), graph = T)
AMAgrowthCurve(mydata.long, treatment = "0.5", ylim = c(0,1), strain = c("EV", "PEN3", "PDR12"),width = 7, height = 5, colours = c("#000000", "#FF3853", "#00BBFF"), graph = T)


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


DoseResponseCurve = function(df, # a data frame in long format
                          time, # the time point at which growth will be measured
                          control = "0", # the reference treatment to which everything will be set at 100%
                          strain, # strains you are interested in
                          colours = NA, # the specified colours, can be a list
                          width = 5, ylim = c(0,NA), height = 5, graph = F){ # graphing information
  if (length(colours) != length(strain)){
    colours = 1:length(strain)
  }
  library(plyr)
  ref = controlSummary(df, varName = "OD600", control = control)
  df = df[df$time == time,]
  df = df[df$strain %in% strain,]
  df$strain = factor(df$strain, levels = strain)
  df = df[grepl("^{0,1}[0-9]{0,}.{0,1}[0-9]{1,}$", df$treatment),] # only numeric treatments as the x-axis will be continuous
  #df$treatment = as.numeric(df$treatment)
  #data[data$strain == "EV", ]$OD600 = (data[data$strain == "EV", ]$OD600 - ref[ref$strain =="EV" &ref$time==0,]$OD600)/(ref[ref$strain == "EV"& ref$time == select_time,]$OD600- ref[ref$strain =="EV" &ref$time==0,]$OD600) * 100
  for(i in levels(as.factor(df$strain))){
    df[df$strain == i, ]$OD600 = (df[df$strain == i, ]$OD600 - ref[ref$strain ==i &ref$time==0,]$OD600)/(ref[ref$strain == i& ref$time == time,]$OD600- ref[ref$strain ==i &ref$time==0,]$OD600) * 100
    
  }
  detach(package:plyr)
  
  df.summary = df %>% group_by(strain, treatment) %>% summarise(growth = mean(OD600), stdev = sd(OD600))
  p = ggplot(data = df.summary, aes(x = treatment, y = growth, group = strain, color = strain)) + geom_line( size = 1)+
    geom_errorbar(aes(ymin=growth-stdev, ymax=growth+stdev), width=.1, size = 1) +
    geom_point(size = 4) +
    xlab("Dose (mM)") + ylab("Relative growth (% of control)") +
    #scale_x_continuous(limits = c(-2,0.5)) +
    scale_y_continuous(limits = ylim, expand = c(0,0)) +
    #scale_fill_manual(values = colours) +
    scale_color_manual(values = colours) +
    # scale_linetype_manual(values = c(1,2,2)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size=1),
          axis.title.x=element_text(size = 20),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length.x = unit(5, "points"),
          axis.ticks.length.y = unit(5,"points") ,
          axis.title.y = element_text(size = 20),
          axis.text = element_text(colour = "black", size=15),
          axis.text.x = element_text(vjust = -0.25),
          legend.position = "none",
          plot.margin = unit(c(20,0,10,0), "points")
    )
  if(graph ==T){
    exptID = readline(prompt = "Enter the experiment ID: ")
    ggsave(file = paste("DoseResponseCurve", paste0(time, "h"), paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}

DoseResponseCurve(pdr.24.2, control = "UN",strain = c("EV", "PEN3", "PDR12"), time ="48", colours = c("#000000", "#FF3853", "#00BBFF"), ylim = c(-5,110),graph = F,height = 5, width = 7)
DoseResponseCurve(mydata.long, control = "0",strain = c("EV", "PEN3", "PDR12"), time ="", colours = c("#000000", "#FF3853", "#00BBFF"), ylim = c(-5,110),graph = T,height = 5, width = 7)

# Basically the same graph as the dose response curve but a bar graph with the points and you can pick 
relGrowthBar = function(df, # a data frame in long format
                             time, # the time point at which growth will be measured
                             control = "0", # the reference treatment to which everything will be set at 100%
                             strain, # strains you are interested in
                             treatment = NA,
                             colours = NA, # the specified colours, can be a list
                             width = 7, ylim = c(0,NA), barLabs = NA, height = 5, graph = F){ # graphing information
  if (length(colours) != length(strain)){
    colours = 1:length(strain)
  }
  library(plyr)
  ref = controlSummary(df, varName = "OD600", control = control)
  df = df[df$time == time,]
  df = df[grepl("^{0,1}[0-9]{0,}.{0,1}[0-9]{1,}$", df$treatment),] # only numeric treatments as the x-axis will be continuous
  if(is.na(treatment[1]) == F){
    df = df[df$treatment %in% treatment,]
    df$treatment = factor(df$treatment, levels = treatment)
  } else{
    df$treatment = factor(df$treatment)
  }
  df = df[df$strain %in% strain, ]
  df$strain = factor(df$strain, levels = strain)

  #df$treatment = as.numeric(df$treatment)
  #data[data$strain == "EV", ]$OD600 = (data[data$strain == "EV", ]$OD600 - ref[ref$strain =="EV" &ref$time==0,]$OD600)/(ref[ref$strain == "EV"& ref$time == select_time,]$OD600- ref[ref$strain =="EV" &ref$time==0,]$OD600) * 100
  for(i in levels(as.factor(df$strain))){
    df[df$strain == i, ]$OD600 = (df[df$strain == i, ]$OD600 - ref[ref$strain ==i &ref$time==0,]$OD600)/(ref[ref$strain == i& ref$time == time,]$OD600- ref[ref$strain ==i &ref$time==0,]$OD600) * 100

  }
  detach(package:plyr)
  
  if(length(levels(df$treatment)) != length(barLabs)){
    barLabs = levels(df$treatment)
  }
  
  p=ggplot(df, aes(x = treatment, y = OD600, fill = strain, group = strain)) + 
    stat_summary(aes(group = strain), colour = "#000000", fun = mean, geom = 'bar', width = 0.6, size = 1, position = position_dodge(width = 0.8))+ 
    stat_summary( aes(y = OD600, group = strain), fun = mean,
                  fun.min = function(x) {
                    pmax(mean(x) - sd(x), 0, na.rm = T)}, 
                  fun.max = function(x) {mean(x) + sd(x)}, 
                  geom = "errorbar", lty =1 , size =0.75, width = 0.25, position = position_dodge(width = 0.8)) +
    geom_jitter( alpha = 0.4, size = 2, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)) +
    xlab("") +   scale_y_continuous(limits = ylim, expand = c(0,0))+
    ylab(lab.y) +  scale_x_discrete(labels = barLabs) +
    scale_fill_manual(values = colours) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black", size=1),
          axis.title.x=element_blank(),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length.x = unit(5, "points"),
          axis.ticks.length.y = unit(5,"points") ,
          axis.title.y = element_text(size = 20),
          axis.text = element_text(colour = "black", size=15),
          axis.text.x = element_text(vjust = -0.25),
          legend.position = "none",
          plot.margin = unit(c(20,0,10,0), "points")
    )
  
  if(graph ==T){
    exptID = readline(prompt = "Enter the experiment ID: ")
    ggsave(file = paste("RelGrowthBar", paste0(time, "h"), paste0(exptID, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
}
relGrowthBar(pdr.24.3, treatment = c("0", "0.1", "0.2", "0.5", "1", "2", "5"), time = "24", strain = c("EV", "PEN3", "PDR12"),  colours = c("#FFFFFF", "#FF3853", "#00BBFF"), ylim = c(0,110),graph = T,height = 5, width = 7)
