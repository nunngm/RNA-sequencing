##Counting DRGs

t0 = y_up #sets the genes of interest for young plants
t12 = m_up #sets the genes of interested for mature plants
t24 = mock_up
# t0 = names(t0[t0<0.01])
# t12 = names(t12[t12<0.01])
# t24 = names(t24[t24<0.01])

t0_down = y_down
t12_down = m_down
t24_down = mock_down

tp = c("t0", "t12", "t24")
allUp = list(names(t0[t0<0.01]),
             names(t12[t12<0.01]),
             names(t24[t24<0.01])
     )
names(allUp) = tp

allDown = list(names(t0_down[t0_down<0.01]),
               names(t12_down[t12_down<0.01]),
               names(t24_down[t24_down<0.01])
     )
names(allDown) = tp

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


# t0_unique = t0[as.integer(! t0 %in% t12) ==1 & as.integer(! t0 %in% t24) == 1]
# t12_unique = t12[as.integer(! t12 %in% t0) ==1 & as.integer(! t12 %in% t24) == 1]
# t24_unique = t24[as.integer(! t24 %in% t0) ==1 & as.integer(! t24 %in% t12) == 1]
# 
# t0_down_unique = t0_down[as.integer(! t0_down %in% t12_down) ==1 & as.integer(! t0_down %in% t24_down) == 1]
# t12_down_unique = t12_down[as.integer(! t12_down %in% t0_down) ==1 & as.integer(! t12_down %in% t24_down) == 1]
# t24_down_unique = t24_down[as.integer(! t24_down %in% t0_down) ==1 & as.integer(! t24_down %in% t12_down) == 1]

## Non-unique genes
# t0_nonunique = t0[as.integer( t0 %in% t12) ==1 | as.integer( t0 %in% t24) == 1]
# t12_nonunique = t12[as.integer( t12 %in% t0) ==1 | as.integer( t12 %in% t24) == 1]
# t24_nonunique = t12[as.integer( t24 %in% t0) ==1 | as.integer( t24 %in% t12) == 1]

#Create totals
mat = matrix("",nrow = 6, ncol = 3)
colnames(mat) = c("Total DEGs", "Earliest Only", "Unique DEGs"
                    #,"Non-Unique DEGs"
     )
tp = c("t0.25", "t12", "t24")
rownames(mat) = c(paste0(tp,"_up"), paste0(tp,"_dwn"))
mat[1:3,1] = lengths(allUp)
mat[1:3,2] = lengths(allUp_earliest)
mat[1:3,3] = lengths(allUp_unique)
mat[4:6,1] = lengths(allDown)
mat[4:6,2] = lengths(allDown_earliest)
mat[4:6,3] = lengths(allDown_unique)

##Output lists
write.csv(mat, "temp.csv")
