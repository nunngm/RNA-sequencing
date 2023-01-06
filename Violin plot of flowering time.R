## Violin plot of flowering time

library(tidyr)
library(ggplot2)
library(dplyr)
library(svglite)
library(agricolae)
library(car)
library(ggthemes)

df = read.table("clipboard", sep = "\t", header=T)
df$genotype = factor(df$genotype, levels = c("Col-0", "npr1-1", "npr4-4D", "n1n4"))
barLabs = c("Col-0", expression(italic("npr1-1")), expression(italic("npr4-4D")), expression(italic("npr1-1\nnpr4-4D")))

p = ggplot(data = df, aes(x = genotype, y = florAge, fill = genotype)) +
  geom_dotplot(binaxis='y',stackdir = "center", dotsize=0.75, binwidth = 1) + 
  geom_violin(trim = F, scale = "width", adjust = 1.2, alpha = 0.4, draw_quantiles = 0.5) +
  scale_x_discrete()+ ylim(0,NA)+
  ylab("Flowering time (dpg)") +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        strip.text.x = element_text(size = 15),
        strip.background = element_rect(colour = "black", fill = "#FFFFFF", size = 1),
        axis.line = element_line(colour = "black", size=0),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank()),
        axis.ticks=element_line(colour = "black", size =1),
        axis.ticks.length = unit(5,"points") ,
        axis.title.y = element_text(size=15),
        axis.text = element_text(color = "black", size=15))
ggsave(filename = "ARR-NPR-22-2_Flower2.svg", plot = p, width = 7, height = 7)
anovaModel = aov(florAge ~ genotype, data = df)
HSD.test(anovaModel, alpha=0.05, "genotype", console=F)$groups

ggplot(., aes(x=hpi:treatment, y=inlier, fill = treatment))
