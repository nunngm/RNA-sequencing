## Graphing co-localisation data
### A note before beginning: Correlation = Pearson's correlation, co-occurance = meander's split co-occurance (where cooccuranceAB is GFP-> stain coocurrance and cooccuranceBA is stain->GFP cooccurance), ICQ = Li's Intensity Correlation Quotient

library(ggplot2)
library(dplyr)
library(agricolae)

mydata= read.table(file= "clipboard",sep= "\t",header =T)
data = mydata

colocGraph = function(data, stain = c("CFW", "ConA"), colocalisation = c("all", "correlation", "cooccuranceAB", "cooccuranceBA", "ICQ"), strain.levels = NA, width = 5, height = 4, print2file = F){
  if (stain == "CFW"||stain == "ConA"){
    df = data[data$remove != "both" & data$remove != stain ,c(rep(T,7), grepl(pattern = stain, x = colnames(data)[8:ncol(data)]))] #removed the samples with bad staining and removes the columns with the other stains data
  } else{
    stop("Set stain to either 'CFW' or 'ConA'")
  }
  colnames(df) = sub(pattern = paste0(stain, "_"), replacement = "" , x = colnames(df)) #renames the columns to remove the stain prefix
  
  if (is.na(strain.levels) == F){
    df$strain = factor(df$strain, levels = strain.levels)
  } else{
    df$strain = as.factor(df$strain)
  }
  df$slide = as.numeric(df$slide)
  df$correlation = as.numeric(df$correlation)
  df$cooccuranceAB = as.numeric(df$cooccuranceAB)
  df$cooccuranceBA = as.numeric(df$cooccuranceBA)
  df$ICQ = as.numeric(df$ICQ)
  
  if(length(colocalisation) >1){
    stop("Please pick a viable option for 'colocalisation'")
  }
  if(colocalisation == "correlation" || colocalisation == "all"){
    anovaModel = aov(correlation ~ strain, data = df)
    print(HSD.test(anovaModel, alpha=0.05, "strain", console=F)$groups)
    
    ggplot(df, aes(x = strain, y = correlation, fill = strain, group = strain)) +
      geom_boxplot(position = position_dodge(width = 0.9), width = 0.8, size = 0.75, outlier.alpha = 0) +
      geom_jitter( alpha = 0.4, size = 2, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8)) + 
      xlab("Strain") +   scale_y_continuous(limits = ylim, expand = c(0,0)) + ylab("Pearson's correlation") +  scale_x_discrete(labels = strain.levels) +
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
            axis.title.y = element_text(size=15),
            axis.text = element_text(colour = "black", size=15),
            axis.text.x = element_text(vjust = -0.25),
            legend.position = "none")
  }
  summary(df$strain)
}