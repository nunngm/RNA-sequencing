## Graphing co-localisation data
### A note before beginning: Correlation = Pearson's correlation, co-occurance = meander's split co-occurance (where cooccuranceAB is GFP-> stain coocurrance and cooccuranceBA is stain->GFP cooccurance), ICQ = Li's Intensity Correlation Quotient

library(ggplot2)
library(dplyr)
library(agricolae)

#laptop directory
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")

mydata= read.table(file= "clipboard",sep= "\t",header =T)
data = mydata

colocGraph = function(data, 
                      stain = c("CFW", "ConA"), 
                      colocalisation = c("correlation", "cooccuranceAB", "cooccuranceBA", "ICQ"), 
                      strain.levels, strain.labels = "", width = 6, height = 7, print2file = F){
  if (stain == "CFW"||stain == "ConA"){
    df = data[data$remove != "both" & data$remove != stain ,c(rep(T,7), grepl(pattern = stain, x = colnames(data)[8:ncol(data)]))] #removed the samples with bad staining and removes the columns with the other stains data
  } else{
    stop("Set stain to either 'CFW' or 'ConA'")
  }
  colnames(df) = sub(pattern = paste0(stain, "_"), replacement = "" , x = colnames(df)) #renames the columns to remove the stain prefix
  
  if (length(strain.labels)==1){
    strain.labels = strain.levels
  }
  
  if (length(strain.levels)>1){
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
  if(colocalisation == "correlation" ){
    y.lab = "Pearson's correlation"
    ylim = c(0,1)
  } else if(colocalisation == "cooccuranceAB") {
    y.lab = "Co-occurance (A->B)"
    ylim = c(0,1)
  } else if( colocalisation == "cooccuranceBA"){
    y.lab = "Co-occurance (B->A)"
    ylim = c(0,1)
  } else if (colocalisation == "ICQ"){
    y.lab = "Li's Intensity Correlation Quotient"
    ylim = c(-0.5,0.5)
  }
    
  anovaModel = aov(get(colocalisation) ~ strain, data = df)
  print(HSD.test(anovaModel, alpha=0.05, "strain", console=F)$groups)
    
  p = ggplot(df, aes(x = strain, y = get(colocalisation), fill = strain, group = strain)) +
      geom_boxplot(position = position_dodge(width = 0.9), width = 0.8, size = 0.75, outlier.alpha = 0) +
      geom_jitter( alpha = 0.4, size = 2, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8)) + 
      xlab("") +   scale_y_continuous(limits = ylim, expand = c(0,0)) + ylab(y.lab) +  scale_x_discrete(labels = slabels) +
      scale_fill_manual(values = c("white", "white", "white", "white")) +
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
            axis.text.x = element_text(vjust = -0.25, face = "italic"),
            legend.position = "none")
  if (print2file ==T){
    exptID = readline(prompt = "Enter experiment ID: ")
    ggsave(file = paste(exptID, stain, paste0(colocalisation, ".svg"), sep = "_"), plot = p, width = width, height = height)
  } else{
    p
  }
  
}
slabels = c("Pst\n", "Pst ΔalgD\n","Pst ΔalgU\nmucAB", "Pst ΔalgU\nmucAB ΔalgD")
slevels = c("Pst", "Pst ΔalgD","Pst ΔalgU mucAB", "Pst ΔalgU mucAB ΔalgD")

colocGraph(data, stain = "ConA", colocalisation = "cooccuranceBA", strain.labels = slabels, strain.levels = slevels, print2file = T)
