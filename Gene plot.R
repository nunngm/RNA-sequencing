library(genemodel)
library(gggenomes)

PCR1 = read.table("clipboard", sep = "\t", header =T)
PCR1[14,1] = "3\' utr"
genemodel.plot(model = PCR1,start = 5132506, bpstop = 5133860, orientation = "reverse", xaxis= T)

genemodel.plot(model = PCR1temp, start = 5132506, bpstop = 5133960, orientation = "reverse")

genemodel.plot(model = PCR1,start =0, bpstop = 20000,orientation = "reverse", xaxis= F)
mutation.plot(5133955,5133920, text = "")

makeGenePlot = function(gene, xlim = NULL, width = 8, height = 1, bpstart = 0, direction = "forward"){
  gene = cbind(gene, t(as.data.frame(strsplit(gene$coordinates, split = "-"))))
  colnames(gene) = c("type", "coordinates", "start", "end")
  if (direction == "forward"){
    gene$start = as.numeric(gene$start) + bpstart
    gene$end = as.numeric(gene$end) + bpstart
  } else if(direction == "reverse"){
    gene$start = bpstart - as.numeric(gene$start)
    gene$end = bpstart - as.numeric(gene$end)
  } else{
    stop("Direction needs to be either 'forward' or 'reverse'")
  }


  d.exons = data.frame(x1 = gene[gene$type == "exon",]$start, x2 = gene[gene$type == "exon",]$end, y1 = 0, y2 = 1)
  introns = gene[gene$type == "intron",]
  introns$start = introns$start-1
  introns$end = introns$end+1
  d.introns.up = data.frame(x1 = introns$start, x2 =(introns$start+ ((introns$end-introns$start)/2)), y1 = 0.5, y2 = 1)
  d.introns.down = data.frame(x1 = (introns$start+ ((introns$end-introns$start)/2)), x2 = introns$end, y1 = 1, y2 = 0.5)
  d.cds = data.frame(x1 = gene[gene$type == "coding_region",]$start, x2 = gene[gene$type == "coding_region",]$end, y1 = 0, y2 = 1)

  p = ggplot() + 
    scale_x_continuous(name="x") + 
    scale_y_continuous(name="y", expand = c(0.1,0.1)) +
    geom_rect(data=d.exons, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill = "blue",color="black", alpha=0.5) +
    geom_rect(data=d.cds, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill = "red",color="black", alpha=1) +
    geom_segment(data = d.introns.up, mapping = aes(x = x1, y = y1, xend = x2, yend = y2)) +
    geom_segment(data = d.introns.down, mapping = aes(x = x1, y = y1, xend = x2, yend = y2)) +
    coord_cartesian( xlim = xlim) +
    geom_vline(xintercept = c(5133955,5133920), colour = "blue", linetype = "longdash")+
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_blank(),
          strip.text.x = element_text(size = 5),
          strip.background = element_rect(colour = "black", fill = "#FFFFFF", size = 1),
          axis.line.x = element_line(colour = "black", size=1),
          axis.line.y = element_blank(),
          axis.title.x=element_blank(),
          #axis.text.x=element_blank()),
          axis.ticks=element_line(colour = "black", size =1),
          axis.ticks.length.x = unit(5,"points") ,
          axis.ticks.length.y = unit(0,"points"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(color = "black", size=10),
          axis.text.y = element_blank())
  
  # if (direction == "reverse"){
  #   p = p+scale_x_reverse()
  # }
  
  geneName = readline(prompt = "enter gene name: ")
  #p
  ggsave(file = paste(geneName, "_geneModel.svg", sep = "_"), plot = p, width = width, height = height)
}


