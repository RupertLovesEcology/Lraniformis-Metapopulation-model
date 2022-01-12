## Remove everything
rm(list = ls())


AllWintSpr <- read.csv("C:/Workspace/18_8_21 Calc Annual Trans.csv", header = T, sep = ",", dec = ".")
AllWintSpr <- AllWintSpr[,-1]
AllWintSpr <- as.numeric(AllWintSpr)



AnnualTotals <- array(data = NA, dim = 83 * 200)
dim(AnnualTotals) <- c(200,83)

Yr = 1924
inc <- 0

for (i in 1:83) {
#  print("starting 1")
  Yr = 1924
  
  Yr <- Yr + i
  cat("year is ", Yr, "\n")
  counter <- 1
  
for (j in 1:nrow(AllWintSpr)) {
# cat("starting 2 j is ", j, "\n")
  if (AllWintSpr[j,1] == Yr) {
#    print("True")
    AnnualTotals[counter,i] <- AllWintSpr[j,2]
    counter <- counter + 1
    cat("counter is ", counter, "\n")
  }
}
#  cat("nrow AnnualTotals for ", Yr, " is ", nrow )
}

Yr <- 1924
 sills <- array(data = NA, dim = 83*2)
 dim(sills) <- c(83,2)
SortTotals <- AnnualTotals 

for (i in 1:83) {
  A <- AnnualTotals[,i]
  SortTotals[,i] <- sort(A, decreasing = T, na.last = T)
}

SortTotals <- SortTotals-0.1


# use the sorted array to determine the threshold  
for (i in 1:83) {
  a <- SortTotals[1,i]
  F <- SortTotals[,i]
  for (f in 1:200) {
  
  if (length(which(as.vector(F) > a)) < 10) {
    a <- a-f/100
  } else {
    break
  }
  cat("a is ", a, "\n")
    }
## Store the 10th highest in row 199
    SortTotals[199,i] <- a
## Store the rounded versions in row 200 for AnnualTransitionsofExisitngSills.R and AnnualTransitionPr.csv   
    SortTotals[200,i] <- round_any(SortTotals[199,i],0.5)
}
 
write.csv(SortTotals,"C:/Workspace/obsAnnualSills.csv",row.names = FALSE)

DroughtAppend <- read.csv("C:/Workspace/DroughtYears.csv", header = T, sep = ",", dec = ".")

for (i in 1:23) {
  DroughtAppend[i,3] <- round_any(DroughtAppend[i,2],0.5)
}
write.csv(DroughtAppend,"C:/Workspace/DroughtYears.csv",row.names = FALSE)
 