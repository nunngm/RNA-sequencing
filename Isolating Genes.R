##Counting DRGs
res = compare.group(age = "y")
tp = c("0.25", "12", "24")
filebase = "y.xlsx"
comp = c("Y.Pst>Y.Mock", "Y.Mock>Y.Pst")
reg = c("Up", "Down")
pval = 0.05
lfc_cutoff = 0

allGenes = lapply(res, function(x){
          up = x[x$log2FoldChange > lfc_cutoff, ]$padj #collect genes where fold change is positive (up-regulated)
          names(up) = rownames(x[x$log2FoldChange > lfc_cutoff, ]) #collect the genenames
          up = up[complete.cases(up)]
          up = names(up[up < pval])
          down = x[x$log2FoldChange < lfc_cutoff, ]$padj #collect genes where fold change is positive (up-regulated)
          names(down) <- rownames(x[x$log2FoldChange < lfc_cutoff, ]) #collect the genenames
          down <- down[complete.cases(down)]
          down = names(down[down < pval])
          up = list(up, down)
          names(up) = c("up", "down")
          return(up)
})

allUp = lapply(allGenes, `[[`, 1)
allDown = lapply(allGenes, `[[`, 2)

#Create earliest expression lists
allUp_earliest = list(
                    allUp[[1]],
                    allUp[[2]][ allUp[[2]] %in% allUp[[1]] == F],
                    allUp[[3]][allUp[[3]] %in% allUp[[1]] == F & allUp[[3]] %in% allUp[[2]] == F]
)

allDown_earliest = list(
                    allDown[[1]],
                    allDown[[2]][ allDown[[2]] %in% allDown[[1]] == F],
                    allDown[[3]][allDown[[3]] %in% allDown[[1]] == F & allDown[[3]] %in% allDown[[2]] == F]
)

#Create unique expression list
allUp_unique = list(
                    allUp[[1]][allUp[[1]] %in% allUp[[2]] == F & allUp[[1]] %in% allUp[[3]] == F],
                    allUp[[2]][ allUp[[2]] %in% allUp[[1]] == F & allUp[[2]] %in% allUp[[3]] == F],
                    allUp[[3]][allUp[[3]] %in% allUp[[1]] == F & allUp[[3]] %in% allUp[[2]] == F]
)

allDown_unique = list(
                    allDown[[1]][allDown[[1]] %in% allDown[[2]] == F & allDown[[1]] %in% allDown[[3]] == F],
                    allDown[[2]][ allDown[[2]] %in% allDown[[1]] == F & allDown[[2]] %in% allDown[[3]] == F],
                    allDown[[3]][allDown[[3]] %in% allDown[[1]] == F & allDown[[3]] %in% allDown[[2]] == F]

)

allUp_nonunique = list(
                    allUp[[1]][allUp[[1]] %in% allUp[[2]] ==T | allUp[[1]] %in% allUp[[3]] ==T],
                    allUp[[2]][allUp[[2]] %in% allUp[[1]] ==T | allUp[[2]] %in% allUp[[3]] ==T],
                    allUp[[3]][allUp[[3]] %in% allUp[[1]] ==T | allUp[[3]] %in% allUp[[2]] ==T]
)

allDown_nonunique = list(
                    allDown[[1]][allDown[[1]] %in% allDown[[2]] ==T | allDown[[1]] %in% allDown[[3]] ==T],
                    allDown[[2]][allDown[[2]] %in% allDown[[1]] ==T | allDown[[2]] %in% allDown[[3]] ==T],
                    allDown[[3]][allDown[[3]] %in% allDown[[1]] ==T | allDown[[3]] %in% allDown[[2]] ==T]
)

#Variables


#Create totals
mat = matrix("",nrow = 6, ncol = 4)
colnames(mat) = c("Total DEGs", "Earliest Only", "Unique DEGs"
                    ,"Non-Unique DEGs"
     )

rownames(mat) = c(paste0(tp,"_up"), paste0(tp,"_dwn"))
mat[1:3, 1] = lengths(allUp)
mat[1:3, 2] = lengths(allUp_earliest)
mat[1:3, 3] = lengths(allUp_unique)
mat[1:3, 4] = lengths(allUp_nonunique)
mat[4:6, 1] = lengths(allDown)
mat[4:6, 2] = lengths(allDown_earliest)
mat[4:6, 3] = lengths(allDown_unique)
mat[4:6, 4] = lengths(allDown_nonunique)

##Output lists

write.xlsx(mat, filebase, sheetName = "Counts", col.names = T, row.names = T)
for(j in 1:2){
     write.xlsx(0,paste0(reg[j], "_", filebase))
     for(i in 1:3){
          time = strsplit(tp, ".", fixed = T)[[i]][1]
          
          df = genes.info(get(paste0("all", reg[j]))[[i]], time)
          write.xlsx(df , paste0(reg[j], "_", filebase), sheetName = paste0(comp[j], " ", tp[i]," hpi") ,col.names = T, row.names = F, append = T)
          df = genes.info(get(paste0("all", reg[j], "_earliest"))[[i]], time)
          write.xlsx(df , paste0(reg[j], "_", filebase), sheetName = paste0(comp[j], " earliest ", tp[i]," hpi") ,col.names = T, row.names = F, append = T)
          df = genes.info(get(paste0("all", reg[j], "_unique"))[[i]], time)
          write.xlsx(df , paste0(reg[j], "_", filebase), sheetName = paste0(comp[j], " unique ", tp[i]," hpi") ,col.names = T, row.names = F, append = T)
          df = genes.info(get(paste0("all", reg[j], "_nonunique"))[[i]], time)
          write.xlsx(df , paste0(reg[j], "_", filebase), sheetName = paste0(comp[j], " nonunique ", tp[i]," hpi") ,col.names = T, row.names = F, append = T)
     }
     removeSheet(paste0(reg[j], "_", filebase))
}



#______________________Go no further__________________________________________

df = genes.info(allUp[[1]], "0")
write.xlsx(df , filebase, sheetName = paste0(comp ," 0.25 hpi") ,col.names = T, row.names = F)
df = genes.info(allUp[[2]], "12")
write.xlsx(df , filebase, sheetName = paste0(comp, " 12 hpi") ,col.names = T, row.names = F, append = T)
df = genes.info(allUp[[3]], "24")
write.xlsx(df , filebase, sheetName = paste0(comp, " 24 hpi") ,col.names = T, row.names = F, append = T)

df = genes.info(allUp_earliest[[1]], "0")
write.xlsx(df , filebase, sheetName = paste0(comp, " early 0.25 hpi"),col.names = T, row.names = F, append = T)
df = genes.info(allUp_earliest[[2]], "12")
write.xlsx(df , filebase, sheetName = paste0(comp, " early 12 hpi"),col.names = T, row.names = F, append = T)
df = genes.info(allUp_earliest[[3]], "24")
write.xlsx(df , filebase, sheetName = paste0(comp, " early 24 hpi"),col.names = T, row.names = F, append = T)



write.xlsx(allUp[[2]],"Young plants.xlsx",sheetName = "Y.Pst>Y>Mock 12 hpi",col.names = F, row.names = F,append = T)
write.xlsx(allUp[[3]],"Young plants.xlsx",sheetName = "Y.Pst>Y>Mock 24 hpi",col.names = F, row.names = F,append = T)
write.xlsx(allUp_earliest[[1]],"Young plants.xlsx",sheetName = "Unique 0.25 hpi",col.names = F, row.names = F,append = T)