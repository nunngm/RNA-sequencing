

library(tidyr)
library(ggplot2)
library(dplyr)
library(svglite)
library(agricolae)
library(car)
library(ggthemes)
library(ggbeeswarm)

pal = colorRampPalette(c( "white", "#54B031"))

df = read.table("clipboard", sep = "\t", header=T)
df$genotype = factor(df$genotype, levels = c("Col-0", "npr1-1", "npr4-4D", "n1n4"))

RD.graph = function(data, colours = c(1,2,3,4), exptID = "temp", width =5, height = 4, graph = F){
  anovaModel = aov(diameter ~ genotype, data = data)
  print(HSD.test(anovaModel, alpha=0.05, "genotype", console=F)$groups)
  
  p = ggplot(data = data, aes(x = genotype, y = diameter, fill = genotype)) +
    #geom_jitter(aes(colour = genotype),  position = position_jitterdodge(jitter.width = 1)) +
    #geom_dotplot(binaxis='y',stackdir = "center", dotsize = 1, binwidth = 0.1) + 
    #geom_beeswarm() +
    geom_quasirandom(aes(colour = genotype), method = "quasirandom") +
    geom_violin(aes(fill = genotype),trim = F, scale = "width", adjust = 1.2, alpha = 0.4, draw_quantiles = 0.5) +
    scale_x_discrete()+ ylim(0,NA)+
    ylab("Rosette diameter(cm)")+ theme(panel.grid.major = element_blank(),  
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 legend.position = "none",
                 panel.border = element_rect(colour = "black", fill = NA, size = 1),
                 strip.text.x = element_text(size = 15),
                 strip.background = element_rect(colour = "black", fill = "#FFFFFF", size = 1),
                 axis.line = element_line(colour = "black", size=0),
                 axis.title.x=element_blank(),
                 axis.text.x=element_text(face = c("plain", "italic", "italic", "italic")),
                 axis.ticks=element_line(colour = "black", size =1),
                 axis.ticks.length = unit(5,"points") ,
                 axis.title.y = element_text(size=15),
                 axis.text = element_text(color = "black", size=15)) + scale_colour_manual(values = colours) +scale_fill_manual(values = colours)
  if (graph == T){
    svg(filename = paste0("RD_", exptID, ".svg"), width = width, height = height )
    p
    dev.off()
    #ggsave(filename = paste0("RD_", exptID, ".svg"),width = width, height = height)
  } else{
    p
  }
}

ggsave(filename = "ARR-NPR-22-2_growth.svg", plot = p, width = 10, height = 8)

