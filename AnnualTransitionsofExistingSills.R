


source("C:/workspace/math0286/R/win-library/3.6/matrixOperators.r")
setwd("C:/workspace/math0286/R")
.libPaths("C:/workspace/math0286/R/win-library/4")




## Remove everything
rm(list = ls())

##Libraries
library(plyr)


StateTrans <- read.csv("C:/Workspace/AnnualTransitionPrV2.csv", header = T, sep = ",", dec = ".")
sillForecast <- read.csv("C:/Workspace/MetapopHydrologyForecast.csv", header = TRUE, sep = ",", dec = ".")
#sillForecast <- sillForecast[,2]
sillForecast <- sillForecast - 0.1
## NOTE *************************************** 
## SILLFORESCAST HAS HAD THE 0.1 DEDUCTED ABOVE

max <- length(sillForecast)

Set <- sillForecast[1:10,1]
Next <- 21

setsize <- 60

## This code groups the sillForecasts into meaningful sets of 60 (keeping inmind we have alreday done 15 burn in and 23 drought) 
## This means we can have 10000 sets to work with

Testi <- array(data = NA, dim = 244)
dim(Testi) <- c(61,4)

Testi[1,1] <- 9
Testi[1,2] <- 8.949
Testi[1,3] <- 8.949
Testi[1,4] <- 1377.024



#Making 10000 sequences of 60 years
HydroSets <- array(data=NA, dim = 60*10000)
dim(HydroSets) <- c(10000,60)

#Making 10000 sequences of 60 years of the annual wetness (every winter and spring wetness )
WetSets <- array(data=NA, dim = 60*10000)
dim(WetSets) <- c(10000,60)



## PLay code for transitions

#startval <- round_any(startval,0.5)


#  MUST BE ONE LOOP BECAUSE IT AFFECTS ITS OWN TRANSITIONS
# Step 1: Make the coarse sequence 


rowCount <- 1
subSet <- array(data = NA, dim = 20)
dim(subSet) <- c(10,2)
subSet[1:10,1] <- sillForecast[1:10,2]
subSet[1:10,2] <- sillForecast[1:10,5]
counter <- 11

for (iter in 1:10000) {
for (i in 1:60) {
  
  
  
#cat("START Testi[i,] is ", Testi[i,], "\n")
#  cat("START Testi[i+1,] is ", Testi[i+1,], "\n")
 
  
  if (Testi[i,1] < 6) { Testi[i,1] <- 6 }
  if (Testi[i,1] > 16) { Testi[i,1] <- 16 }
 trans <- runif(1,0,1)
 rowref <- match(c(Testi[i,1]),StateTrans[,1])
 if (trans <= StateTrans[rowref,2]) {
   Testi[i + 1,2] <- runif(1,(Testi[i,3] - StateTrans[rowref,5]),Testi[i,3])
 } else if (trans > StateTrans[rowref,2] && trans <= StateTrans[rowref,3])  {
   Testi[i + 1,2] <- runif(1,Testi[i,3],(Testi[i,3] + StateTrans[rowref,6]))
 } else if (trans > StateTrans[rowref,3] && trans <= StateTrans[rowref,4]) {
   Testi[i + 1,2] <- Testi[i,3]
 } else {
   stop("shit fcked")
 }
  #find the index closest to the predicted value from the first 10 forecasted values
  ind <- which.min(abs(subSet-Testi[i+1,2]))
  Testi[i+1,3] <- subSet[ind,1]
  Testi[i+1,4] <- subSet[ind,2]
  
  # store it in column 3 of Testi then replace the one we used with the next forecasted value
  if (counter <= nrow(sillForecast)) {
  subSet[ind,1] <- sillForecast[counter,2]
  subSet[ind,2] <- sillForecast[counter,5]
  } else {
    subSet[ind,1] <- NA
    subSet[ind,2] <- NA
  }
  counter <- counter + 1
  Testi[i+1,1] <- round_any(Testi[i,3],0.5)
#  cat("END Testi[i,] is ", Testi[i,], "\n")
#  cat("END Testi[i+1,] is ", Testi[i+1,], "\n")
  HydroSets[rowCount,i] <- Testi[i+1,3]
  if (i >= 2) {
  
    WetSets[rowCount,i-1] <- Testi[i+1,4]  
    }
 
  if (i == 60) { 
    rowCount <- rowCount + 1
    print(rowCount)
    }
  }
}


write.csv(HydroSets,"C:/Workspace/OrderedHydroForecast.csv",row.names = FALSE)
write.csv(WetSets,"C:/Workspace/OrderedAnnualWetness.csv",row.names = FALSE)


