## Absolute quantification
### Based on Rutledge and Stewart, 2008
# Optical calibration factor (FU/ng dsDNA) for our kit (Biorad Universal SYBR green) = 957 FU/ng dsDNA
# F0 = calculated initial fluorescence, M0 = calculated initial mass of target double-stranded cDNA, N0 number of initial transcripts, OCF = optical calibration factor, As = amplicon size
# The N0 is scaled by amplicon size irrespective of GC content
library (qpcR)
mydata= read.table(file= "clipboard",sep= "\t",header =T)
samples = as.character(mydata[1,2:ncol(mydata)])
mydata = mydata[2:nrow(mydata),]
groups = t(as.data.frame(strsplit(samples, split = "_")))
colnames(groups) = 
#laptop directory
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-qRT-PCR\\Primer Efficiency test")

#Desktop directory
setwd("C:\\Users\\garrett\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-qRT-PCR\\Primer Efficiency test")

curve = pcrfit(mydata, fluo= 2, cyc = 1, model = l4)
plot(curve)
groups = 


E = efficiency(curves, type = "cpD2",thresh = Ft)
Ct = E$cpT
Emax = E$eff/2
Ft = 187
F0 = Ft/(Emax+1)^Ct
OCF = 957
M0 = F0/OCF
N0 = M0 * 9.1 * 10^11 / As # where 9.1 x 10^11 is the number of base pairs per nanogram of dsDNA.

Ft = 377 # Average fluorescence
OCF = 356 # OCF of a similar 
As = rep(121, times = ncol(mydata)-1)

curves = lapply(2:ncol(mydata), function(x){pcrfit(mydata, cyc = 1, fluo = x)})
E = data.frame(t(sapply(1:(ncol(mydata)-1), simplify = "array", function(x){efficiency(curves[[x]], type = "cpD2", thresh = Ft)})))
Ct = unlist(E$cpT)
Emax = unlist(E$eff)/2
F0 = Ft/(Emax+1)^Ct
M0 = F0/OCF
N0 = M0 * 9.1 * 10^11 / As
N0
