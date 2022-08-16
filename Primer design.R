reverseComplement = function(input, rev = T){
  input = chartr("ATGCatgc","TACGtacg", input)
  if (rev==T){
    reversed = strsplit(as.character(input), split ="")
    input = reversed[[1]][nchar(input):1]
  }
  return(paste(input,collapse = ""))
}


primerFun = function(input){
  input = chartr("ATGC","atgc", input)
  input = unlist(strsplit(as.character(input), split = ""))
  count = 0
  for (i in 1:length(input)){
    if (input[i]== "g" | input[i] == "c"){
      count = count+1
    }
  }
  
  
  print(paste("length:", length(input)))
  print(paste("GC content:", round(count/length(input)*100, digits = 2), "%"))
  print(paste("Melting temp (deg C):", 64.9+41*(count-16.4)/length(input)))
  print(paste("Reverse complement:", chartr("atgc","tacg", paste(input[length(input):1], collapse = ""))))
}
forward ="accttcggaacctttcagcgcacgacaatc"
primerFun(forward)
tbRevd = "CAACACAATCCACCTCGACC"


primerFun(forward)
primerFun(tbRevd)
reverseComplement("GGGTTCACAACCAGGCAC",T)


#----------------------------
# Primer efficiency
library(qpcR)
mydata= read.table(file= "clipboard",sep= "\t",header =T)

#laptop directory
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-qRT-PCR\\Primer Efficiency test")

#Desktop directory
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-qRT-PCR\\Primer Efficiency test")


curve = pcrfit(mydata, cyc = 1, fluo = 2:4, model = l4)
Cy0(curve, plot = T)
plot(curve)

curves = lapply(2:13, function(x){pcrfit(mydata, cyc = 1, fluo = x)})

par(mfrow = c(1,2))
efficiency(curve, type = "cpD2")

sliwin(curves[[12]])


