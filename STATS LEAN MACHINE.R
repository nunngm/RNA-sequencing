#Garrett's Awesome R stat playing thing plus rough graph at the end

#install packages, only needs to be done once per computer
install.packages("car")
install.packages("agricolae") 

#clear the global environment
rm(list=ls()) 


#Reset and read the data
rm(list=ls())
mydata = read.csv(file = file.choose(),header=T) #pick a file
mydata= read.table(file= "clipboard",sep= "\t",header =F) #read from clipboard
mydata <- read.table("AMA-18-2.csv",sep = ",", header=TRUE,row.names =1, check.names = F) #tell R to read the data file and add it to the global environment 
setwd("C:/Users/garre/OneDrive/Documents/Uni Year 4 Cameron Lab/Data/R WD")

colnames(mydata) = c("Y.Col", "Y.ltpg5", "Y.LTPG5OX", "Y.mctp9", "Y.med16", "Y.stmp10", "Y.ugt72e2")

colnames(mydata) = c("Col.un", "col.12", "col.24", "col.48", "ald.un", "ald.12", "ald.24", "ald.48", "fmo.un", "fmo.12", "fmo.24", "fmo.48")

colnames(mydata) = c("Y.Col", "M.Col", "Y.ald1", "M.ald1", "Y.fmo1", "M.fmo1")

colnames(mydata) = c("Y.Col", "Y.pen3-4", "Y.pdr12-3", "Y.p3p12", "M.Col", "M.pen3-4", "M.pdr12-3", "M.p3p12")
colnames(mydata) = c("Y.Col", "Y.pen3-4", "Y.pdr12-3", "M.Col", "M.pen3-4", "M.pdr12-3")

##Set-up for Stats
library(agricolae)
library(car) 
setup= function(md,forANOVA=F){
  numcol= ncol(md)
  mydp=data.frame()
  orgData= data.frame()
  orgData.names = NULL
  mean_cfu=NULL
  numrep=NULL
  for (i in 1:numcol){
    md[,i]= as.numeric(md[,i])
    mydp = rbind(mydp,sd(md[,i],na.rm=T))
    mean_cfu = rbind(mean_cfu, mean(md[,i],na.rm=T))
    x=0
    for (j in 1:nrow(md)){
      orgData = rbind(orgData, md[j,i])
      orgData.names = rbind(orgData.names, names(md)[i])
      x=x+(1-as.integer(is.na(md[j,i])))
    }
    numrep = c(numrep,rep(x,each = j))
    }
  results = NULL
  orgData = cbind(orgData.names,orgData,numrep)
  names(orgData)= c("type","cfu","reps")
  orgData = subset(orgData, reps >1 )
  orgData = delete.na(orgData,"greater")
  rownames(orgData) = 1:nrow(orgData)
  
  mydp = cbind(names(md),mean_cfu,mydp)
  names(mydp)= c("type","cfu","stdev")
  if (forANOVA==T){
    return(orgData)
  }
  results$orgData = orgData
  results$groupData = mydp
  return(results)
}
analysis=setup(mydata,forANOVA=T)

#Automated
AutomatedStats = function(md, noTukey = F){ 
  results = matrix(NA,ncol =4,nrow =2)
  colnames(results) = c("ANOVA","Shapiro-Wilk","Levene","Bartlett")
  rownames(results) = c("P-value","Result")
  anovaModel = aov(cfu ~ type, data=md) #the numeric variable must be first in this expression 
  r = summary(anovaModel)
  results[1,1] = r[[1]]$`Pr(>F)`[1]
  md.lm = lm(cfu ~ type, data = md)
  md.res = resid(md.lm)
  results[1,2]=shapiro.test(md.res)[2]$p.value
  library(car) 
  results[1,3]=leveneTest(cfu ~ type, data=md)[3][1,1]
  results[1,4] = bartlett.test(cfu ~ type, data=md)[3]$p.value
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
  y$groups = HSD.test(anovaModel, alpha=0.05, "type", console=T)$groups#gives letter codes but no p values
  y$stats=results
  return(y)
}
stats = AutomatedStats(analysis)
View(stats$stats)+ View(stats$groups)
write.table(stats$stats,file= "Stats.txt",sep= "\t",col.names =T, row.names = T)
write.table(stats$groups,file= "Groups.txt",sep= "\t",col.names =T, row.names = T)
##Failure, so you must go back and log data
analysis[,2]= log10(analysis[,2]) #Now redo anova and assumptions

##Slowly do ANOVA and assumptions
{ 
  anovaModel = aov(cfu ~ type, data=orgData) #the numeric variable must be first in this expression 
  summary(anovaModel) #show results, NB, before interpreting the results you need to check that your data meet the assumptions of ANOVA, however, to do this you must first run the ANOVA
  #if pr<0.05 this means that there is a difference in the some of the means
  ##Shapiro, pass if p>0.05, null hypothesis is that there is normality
  #not using ANOVA + other things in ".lm"
  orgData.lm = lm(cfu ~ type, data = orgData)
  orgData.res = resid(orgData.lm)
  names(orgData.res)= orgData[is.na(orgData[,2])==F,1]
  shapiro.test(orgData.res)
  rm(orgData.res)
  #using anova
  shapiro.test(resid(anovaModel))
  #tests for homegeneity of variance, null hypothesis: homogeneity of variance 
  library(car) 
  x=leveneTest(orgData$cfu ~ orgData$type, data=orgData)[3][1,1] #passes if p>0.05 thus variance if homogenous
  bartlett.test(cfu ~ type, data=orgData) #if you can run this run this but is more fickle
  #POST-HOC TESTING, If your data passed assumptions and ANOVA
  library(agricolae)
  TukeyHSD(anovaModel) #gives p values for pair-wise comparisons 
  HSD.test(anovaModel, alpha=0.05, "type", console =T) #gives letter codes but no p values
  
}
#reorder tukey

reorderTukey<-function(inV){
  collapsed <- paste(inV,sep="",collapse = "")
  u <- unique(strsplit(collapsed,"")[[1]])
  if(length(u)<2){
    return(inV)
  }
  u <- u[order(u)]
  m <- matrix(nrow=NROW(inV),ncol=length(u))
  m[]<-F
  for(i in 1:length(inV)){
    s <- strsplit(inV[i],"")[[1]]
    index <- match(s,u)
    m[i,index] <- T
  }
  for(i in 1:(length(u)-1)){
    firstColT <- match(T,m[,i])[1] #first row with true in current column
    firstT <- match(T,rowSums(m[,i:length(u)] > 0))[1] #first row with true in rest
    if(firstT < firstColT){
      colT <- match(T,m[firstT,i:length(u)])[1]
      colT <- colT + i - 1 #correct index for leftout columns in match
      tmp <- m[,colT]
      m[,colT] <- m[,i]
      m[,i] <- tmp
    }
  }
  res <- vector(mode = "character", length=length(trt))
  for(i in 1:length(inV)){
    l <- u[m[i,]]
    res[i] <- paste(l,sep="",collapse = "")
  }
  return(res)
}
stats$groups <- stats$groups[rev(rownames(stats$groups)),] #order the result the way you want
stats$M <- reorder(as.character(stats$M))
library(ggplot2)

##young and old plants (make sure to copy formatted based on file "ARR Expr.pzfx":ARR-17-3)
{
rm(list=ls()) 
  
mydata= read.table(file= "clipboard",sep= "\t",header =F,stringsAsFactors = F) #read from clipboard
age=t(mydata[1,2:ncol(mydata)])
rownames(mydata)=mydata[,1]
mydata = t(mydata[2:nrow(mydata),2:ncol(mydata)])
mydata <- as.data.frame(apply(mydata,FUN= as.numeric,MARGIN=c(1,2)))

# stuff = ordered(mydp,levels = c("young","mature"))
# orgData = as.numeric(mydata[1:nrow(mydata),2:ncol(mydata)])
mydp= NULL
orgData= data.frame()
mean_cfu=NULL
mydp.age = NULL

for (i in 1:ncol(mydata)){
  mean_cfu = rbind(mean_cfu,mean(mydata[1:3,i],na.rm=T))
  mean_cfu = rbind(mean_cfu,mean(mydata[4:6,i],na.rm=T))
  mydp = rbind(mydp,sd(mydata[1:3,i],na.rm=T))
  mydp = rbind(mydp,sd(mydata[4:6,i],na.rm =T))
  mydp.age =rbind(mydp.age,age[1])
  mydp.age = rbind(mydp.age,age[4])
}
mydp = cbind(rep(colnames(mydata[,1:ncol(mydata)]),each =2),mean_cfu,mydp,mydp.age)
colnames(mydp) = c("type","cfu","stdev","age")
orgData.names = as.data.frame(rep(colnames(mydata),each = 6))
orgData.age = as.data.frame(rep(age, times = ncol(mydata)))
orgData = unlist(mydata)
names(orgData)= rep_len(NA, length(orgData))
orgData = cbind(orgData.names,orgData,orgData.age)
colnames(orgData) = c("type","cfu","age")

rm(orgData.names,age,orgData.age,mydp.age,i,mean_cfu)
# orgData = as.data.frame(orgData,stringsAsFactor =F)
# orgData[,3]=as.list.factor(orgData[,3])
# orgData[,3] = factor(orgData[,3],levels = c("2","1"), ordered= T, labels = c("young","mature"))
# orgData[,1] = as.factor(orgData[,1])
mixed <- with(orgData, interaction(type, age))
}

#---Before the next step make sure to remove any NA rows from "orgData" the data or it will screw with the ANOVA---#

##ANOVA for young and old
{
  anovaModel = aov(cfu ~ mixed, data = orgData)
  summary(anovaModel)
  #test for normality of the residuals
  orgData.lm = lm(cfu ~ mixed, data = orgData)
  orgData.res = resid(orgData.lm)
  names(orgData.res)= orgData[,1]
  shapiro.test(orgData.res)
  
  shapiro.test(resid(anovaModel)) 
  #pass if p>0.05

  #tests for homegeneity of variance, null hypothesis: homogeneity of variance 
  library(car) 
  leveneTest(cfu ~ mixed, data=orgData) #passes if p>0.05 thus variance if homogenous
  
  bartlett.test(cfu ~ mixed, data=orgData) #if you can run this run this but is more fickle
  
  #POST-HOC TESTING 
  #only if your data satisfy the assumptions and the ANOVA returned a significant result
  library(agricolae)
  TukeyHSD(anovaModel) #gives p values for pair-wise comparisons 
  print(HSD.test(anovaModel, alpha = 0.05,trt ="mixed")) #gives letter codes but no p values
  anovaModel$model$mixed = make.unique(as.character(anovaModel$model$mixed))
  hsd
} 
library(dplyr)
##if it fails
  orgData[,2]=log10(orgData[,2])
##Residual Plot
par(mfrow = c(2,2))
plot(orgData.lm)
par(mfrow = c(1,1))

#obtain some basic outputs  
print(mydata) 
plot(mydata$type,mydata$cfu) 
tapply(mydata$cfu,mydata$type,mean)
tapply(mydata$cfu,mydata$type,sd)
tapply(mydata$cfu,mydata$type,length) 