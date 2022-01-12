############   Trial 4 working with timeseries data to produce forecasts
###  Divide daily height to winter and spring.
## round each daily river height to 10cm increments
## create a discrete markov chain for each value (upper limit always decreases, lower always increases)
## create a StDev for each change to increment
## recreate a new sample for each year then gather the 10th highest
## 
############

source("C:/workspace/math0286/R/win-library/3.6/matrixOperators.r")
setwd("C:/workspace/math0286/R")
.libPaths("C:/workspace/math0286/R/win-library/4")

library(plyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(moments)
library(sn)
library(e1071)

dailyHeight <- read.csv("C:/Workspace/CompHeightL3.csv", header = TRUE, sep = ",", dec = ".")
colnames(dailyHeight)[1] <- 'Day'
colnames(dailyHeight)[4] <- 'Height'
dailyHeight <- dailyHeight[,1:4]
dailyHeight <- ts(dailyHeight)


winter <- array(data = 50, dim = 4 * 1)
dim(winter) <- c(1,4)
colz <- c("WinDay","WinMonth","WinYear","WinHeight")
colnames(winter) <- colz

spring <- array(data = 50, dim = 4 * 1)
dim(spring) <- c(1,4)
colz <- c("SprDay","SprMonth","SprYear","SprHeight")
colnames(spring) <- colz

for (i in 1:nrow(dailyHeight)) {
  if (6<=dailyHeight[i,2] && dailyHeight[i,2] <= 8) {
    winter <- rbind(winter,dailyHeight[i,])
  } else if (9<=dailyHeight[i,2] && dailyHeight[i,2] <= 11) {
    spring <- rbind(spring,dailyHeight[i,])
  }
}

## Construct our calculating columns
## Note for column 6: 1-decrease,2-same,3-increase
#rounded <- round(winter[,4], digits = 1)
winter <- cbind(winter,0)
winter[,5] <-round(winter[,4], digits = 1)

#winter <- cbind(winter,rounded)
winter <- cbind(winter,NA,NA,NA)
colnames(winter)[6] <- 'markovChange'
colnames(winter)[7] <- 'winHeightChngUp'
colnames(winter)[8] <- 'winHeightChngDwn'

rounded <- round(spring[,4], digits = 1)
spring <- cbind(spring,rounded)
spring <- cbind(spring,NA,NA,NA)
colnames(spring)[6] <- 'markovChange'
colnames(spring)[7] <- 'sprHeightChngUp'
colnames(spring)[8] <- 'sprHeightChngDwn'

## Populate our Columns
## NOTE THERE IS CURRENTLY ISSUES HERE: THE 5.4M WINTER RECORDS 1/3 WAS A DECREASE
## i AM CURRENTLY CALCULATING THE CURRENT POINT IN TIME AS A DECREASE OR AN INCREASE NOT THE DESTINATION ??!!

for (i in 2:(nrow(winter)-1)) {
  if (winter[i,3] == winter[i+1,3] && winter[i,5] < winter[i+1,5] ) {
    winter[i,6] <- 3
    winter[i,7] <- (winter[i+1,5] - winter[i,5])
  } else if (winter[i,3] == winter[i+1,3] && winter[i,5] > winter[i+1,5]) {
    winter[i,6] <- 1
    winter[i,8] <- (winter[i,5] - winter[i+1,5])
  } else if (winter[i,3] == winter[i+1,3] && winter[i,5] == winter[i+1,5]) {
    winter[i,6] <- 2
  } else if (winter[i,3] != winter[i+1,3]) {
    winter[i,6] <- NA
  } else {
    stop("this loop has encountered a logical exception you silly chap")
  }
}

for (i in 2:(nrow(spring)-1)) {
  if (spring[i,3] == spring[i+1,3] && spring[i,5] < spring[i+1,5] ) {
    spring[i,6] <- 3
    spring[i,7] <- (spring[i+1,5] - spring[i,5])
  } else if (spring[i,3] == spring[i+1,3] && spring[i,5] > spring[i+1,5]) {
    spring[i,6] <- 1
    spring[i,8] <- (spring[i,5] - spring[i+1,5])
  } else if (spring[i,3] == spring[i+1,3] && spring[i,5] == spring[i+1,5]) {
    spring[i,6] <- 2
  } else if (spring[i,3] != spring[i+1,3]) {
    spring[i,6] <- NA
  } else {
    stop("this loop has encountered a logical exception you silly chap")
  }
}


# Aight lets Start by creating 2 arrays, one for each of Spring and Winter markov chains

min(winter[,5])
max(winter[2:nrow(winter),5])
unique(winter[2:nrow(winter),5])

markovWinter <- array(data = NA, dim = 107*8)
dim(markovWinter) <- c(107,8)
colz <- c("winHeight","PrDown","PrSame","PrUp","DecMin","DecMax","IncMin","IncMax")
colnames(markovWinter) <- colz
markovWinter[1,1] <- 5.4
for (i in 2:107) {
  markovWinter[i,1] <- markovWinter[i-1,1] + 0.1
}
markovWinter[,5] <- 0.1
markovWinter[,7] <- 0.1

min(spring[,5])
max(spring[2:nrow(spring),5])
markovSpring <- array(data = NA, dim = 102*8)
dim(markovSpring) <- c(102,8)
colz <- c("sprHeight","PrDown","PrSame","PrUp","DecMin","DecMax","IncMin","IncMax")
colnames(markovSpring) <- colz
markovSpring[1,1] <- 5.9
for (i in 2:102) {
  markovSpring[i,1] <- markovSpring[i-1,1] + 0.1
}
markovSpring[,5] <- 0.1
markovSpring[,7] <- 0.1



## Make a list of the available heights (winterCalcParams) 
wintCalcParam <- min(winter[,5])
top <- max(winter[2:nrow(winter),5])

## 
 while (max(wintCalcParam) < top) {
wintCalcParam <- append(wintCalcParam,(max(wintCalcParam) + 0.1))
 }

wintCalcParam <- round(wintCalcParam, digits = 1)


## Make calc arrays
calcArray <- array(data = 0, dim = 3)
dim(calcArray) <- c(1,3)
colnames(calcArray) <- c("down","same","up")
inc <- 0.1
dec <- 0.1

## Repeat for Spring
## Make a list of the available heights (sprCalcParam)
sprCalcParam <- min(spring[,5])
top <- max(spring[2:nrow(spring),5])

## 
while (max(sprCalcParam) < top) {
  sprCalcParam <- append(sprCalcParam,(max(sprCalcParam) + 0.1))
}

sprCalcParam <- round(sprCalcParam, digits = 1)











## Populate the markovWinter array using wintCalcParam and winter then save as csv for later use


for (i in 1:length(wintCalcParam)){
 counterDecrease <- 0
#   cat(" I am processing ", wintCalcParam[i], "\n")
   for (j in 2:nrow(winter)) {
          if (winter[j,5] == wintCalcParam[i]) {
       calcArray[1,winter[j,6]] <- calcArray[1,winter[j,6]] + 1
      if (!is.na(winter[j,6]) && winter[j,6] == 1) {
  #     cat("I have found a decrease on row ", i, " sill height ", winter[j,6])
       dec <- append(dec,winter[j,8])
       counterDecrease <- counterDecrease + 1
     }
     if (winter[j,6] == 3 && !is.na(winter[j,6])) {
       inc <- append(inc,winter[j,7])
     }
       } 
 
   }
  cat("populating the matrix for run ", i, "\n")
  cat(" counterDecrease is ", counterDecrease, "\n")
  cat("dec is ", length(dec), " long", "\n")
        ##  populate  markovWinter for saving
   markovWinter[i,2] <- calcArray[1,1]/sum(calcArray[1,1:3])
        markovWinter[i,4] <- calcArray[1,3]/sum(calcArray[1,1:3])
        ## give any bonus to neutral
        markovWinter[i,3] <- 1 - (markovWinter[i,2] + markovWinter[i,4])
        markovWinter[i,6] <- max(dec)
        markovWinter[i,8] <- max(inc)
  
##  reset our holders
        cat("resetting holders for run ", i, "\n")
  calcArray <- array(data = 0, dim = 3)
  dim(calcArray) <- c(1,3)
  colnames(calcArray) <- c("up","down","same")
  inc <- 0.1
  dec <- 0.1 
  
  counterDecrease <- 0
  cat(" counterDecrease is ", counterDecrease, "\n")
  cat("dec is ", length(dec), " long", "\n")
}

markovWinter <- markovWinter[1:106,]

write.csv(markovWinter,'C:/Workspace/markovWinter.csv')


## Double check the  shared arrays are empty
calcArray <- array(data = 0, dim = 3)
dim(calcArray) <- c(1,3)
colnames(calcArray) <- c("down","same","up")
inc <- 0.1
dec <- 0.1


### Repeat for Spring

for (i in 1:102){
  counterDecrease <- 0
  #   cat(" I am processing ", wintCalcParam[i], "\n")
  for (j in 2:nrow(spring)) {
    if (spring[j,5] == sprCalcParam[i]) {
      calcArray[1,spring[j,6]] <- calcArray[1,spring[j,6]] + 1
      if (!is.na(spring[j,6]) && spring[j,6] == 1) {
        #     cat("I have found a decrease on row ", i, " sill height ", winter[j,6])
        dec <- append(dec,spring[j,8])
        counterDecrease <- counterDecrease + 1
      }
      if (spring[j,6] == 3 && !is.na(spring[j,6])) {
        inc <- append(inc,spring[j,7])
      }
    } 
    
  }
  cat("populating the matrix for run ", i, "\n")
  cat(" counterDecrease is ", counterDecrease, "\n")
  cat("dec is ", length(dec), " long", "\n")
  ##  populate  markovSpring for saving
  markovSpring[i,2] <- calcArray[1,1]/sum(calcArray[1,1:3])
  markovSpring[i,4] <- calcArray[1,3]/sum(calcArray[1,1:3])
  ## give any bonus to neutral
  markovSpring[i,3] <- 1 - (markovSpring[i,2] + markovSpring[i,4])
  markovSpring[i,6] <- max(dec)
  markovSpring[i,8] <- max(inc)
  
  ##  reset our holders
  cat("resetting holders for run ", i, "\n")
  calcArray <- array(data = 0, dim = 3)
  dim(calcArray) <- c(1,3)
  colnames(calcArray) <- c("up","down","same")
  inc <- 0.1
  dec <- 0.1 
  
  counterDecrease <- 0
  cat(" counterDecrease is ", counterDecrease, "\n")
  cat("dec is ", length(dec), " long", "\n")
}

write.csv(markovSpring,'C:/Workspace/markovSpring.csv')





## Then lets determine the starting point for each winter
## This requires a mean and SD of the 1/June from each recorded year
winterStart <- 50
for (i in 1:nrow(winter)) {
  if (winter[i,2] == 6 && winter[i,1] == 1) {
    winterStart <- append(winterStart,winter[i,5])
  }
}

## Randomise a start height from the skewed distribution of observed river heights on the first of July
winterStart <- ts(winterStart[2:length(winterStart)])
startMn <- round(mean(winterStart), digits = 1)
startSD <- round(sd(winterStart), digits = 1)
startSkew <- skewness(winterStart)
startKurt <- kurtosis(winterStart)

## note the below can be exported in a more manageable format include the sn package for cp2dp
# then change params to params <- cp2dp(c(7, 0.9, 1.313998, 4.540583), "ST")
params <- cp2dp(c(startMn, startSD, startSkew, startKurt), "ST")
forecastedYear[1,1] <- round(replicate(1, rst(1, dp = params)),digits = 1)


## I have commented this out for now, has more starting heights than I cpould ever need so don't need to rerun
#startYear <- array(data = NA, dim = 100000000)
#dim(startYear) <- c(100000000,1)
#colnames(startYear) <- c("mAHD on First Day of Winter")

#for (i in 1:100000000) {
#  first <- as.numeric(round(replicate(1, rst(1, dp = params)),digits = 1))
#  startYear[i,1] <- first
#  if (i%%1000==0) {
#  cat("Run ", i, "\n")
#  }
#}

# write.csv(startYear,'C:/Workspace/FirstHeightofWinter.csv')
startYear <- read.csv("C:/Workspace/FirstHeightofWinter.csv", header = TRUE, sep = ",", dec = ".")

library(dplyr)
altWint <- read.csv("C:/Workspace/alteredWinter1.csv", header = TRUE, sep = ",", dec = ".")
altSpr <- read.csv("C:/Workspace/alteredSpring1.csv", header = TRUE, sep = ",", dec = ".")

#### now to calculate a timeseries for the new year

## Create dataframes so that the 'filter' function works 
altmarkWintDF <- data.frame(altWint)
altmarkWintDF <- altmarkWintDF[,2:9]
altmarkWintDF[,1] <- round(altmarkWintDF[,1], digits = 1)
altmarkWintDF[,3] <- altmarkWintDF[,3] + altmarkWintDF[,2]
altmarkWintDF[,4] <- altmarkWintDF[,4] + altmarkWintDF[,3]

altmarkSprDF <- data.frame(altSpr)
altmarkSprDF <- altmarkSprDF[,2:9]
altmarkSprDF[,1] <- round(altmarkSprDF[,1], digits = 1)
altmarkSprDF[,3] <- altmarkSprDF[,3] + altmarkSprDF[,2]
altmarkSprDF[,4] <- altmarkSprDF[,4] + altmarkSprDF[,3]

## for starters let's create an array to hold our new series 
# winter - 92 d, spring 91 d - total 183
forecastedYear <- array(data = NA, dim = 183)
dim(forecastedYear) <- c(183,1)
colnames(forecastedYear) <- c("Height")




forecastedYear[1,1] <- as.numeric(startYear[1,2])
forecastedYear[1,1] <- round(forecastedYear[1,1],digits = 1)

## first comes winter (note starting height has been calculated)
for (i in 1:92) {
  currentH <- as.numeric(forecastedYear[i,1])
  currentH <- round(currentH,digits = 1)
 activeRow <-  filter(altmarkWintDF, winHeight == currentH)
 Prchange <- runif(1,0,1)
 if (Prchange <= activeRow[2]) {
 chng <- runif(1,as.numeric(activeRow[5]),as.numeric(activeRow[6]))
     forecastedYear[i+1,1] <- forecastedYear[i,1] - round(chng, digits = 1)  
 } else if (Prchange > activeRow[2] && Prchange <= activeRow[3]) {
   forecastedYear[i+1,1] <- forecastedYear[i,1]
 } else if (Prchange > activeRow[3] && Prchange <= activeRow[4]) {
   chng <- runif(1,as.numeric(activeRow[7]),as.numeric(activeRow[8]))
   forecastedYear[i+1,1] <- forecastedYear[i,1] + round(chng, digits = 1)  
 } else {
   stop("Shit is fucked")
 }
}

## then spring (an additional 91 days)
for (i in 93:182) {
  currentH <- forecastedYear[i,1]
  if (i == 93 && forecastedYear[i,1] < 5.8) { forecastedYear[i,1] <- 5.8 }
  currentH <- round(currentH,digits = 1)
  activeRow <-  filter(altmarkSprDF, sprHeight == currentH) 
  Prchange <- runif(1,0,1)
  if (Prchange <= activeRow[2]) {
    chng <- runif(1,as.numeric(activeRow[5]),as.numeric(activeRow[6]))
    forecastedYear[i+1,1] <- forecastedYear[i,1] - round(chng, digits = 1) 
  } else if (Prchange > activeRow[2] && Prchange <= activeRow[3]) {
    forecastedYear[i+1,1] <- forecastedYear[i,1]
  } else if (Prchange > activeRow[3] && Prchange <= activeRow[4]) {
    chng <- runif(1,as.numeric(activeRow[7]),as.numeric(activeRow[8]))
    forecastedYear[i+1,1] <- forecastedYear[i,1] + round(chng, digits = 1) 
  } else {
    stop("Shit is fucked")
  }
}
plot(forecastedYear)


library(Rfast)

## REPEAT THE ABOVE PROCESS BUT IN A LOOP
## Create an array to hold the annual Sill thresholds and sum of winter and sum of spring
annualSeries <- array(data = NA, dim = 1800000)
dim(annualSeries) <- c(600000,3)
colnames(annualSeries) <- c("sillHeightFilled", "sumofWinter","sumofSpring")

## Creates 150 rows/sets of 1 century 

for (ll in 1:600000) {

  forecastedYear[1,1] <- as.numeric(startYear[ll,2])
  if (forecastedYear[1,1] < 5.3) { forecastedYear[1,1] <- 5.3 }
  if (forecastedYear[1,1] > 15.9) { forecastedYear[1,1] <- 15.9 }
  forecastedYear[1,1] <- round(forecastedYear[1,1],digits = 1)
  
  ## first comes winter (note starting height has been calculated)
  for (i in 1:92) {
    currentH <- as.numeric(forecastedYear[i,1])
    currentH <- round(currentH,digits = 1)
    activeRow <-  filter(altmarkWintDF, winHeight == currentH)
    Prchange <- runif(1,0,1)
    if (Prchange <= activeRow[2]) {
      chng <- runif(1,as.numeric(activeRow[5]),as.numeric(activeRow[6]))
      forecastedYear[i+1,1] <- forecastedYear[i,1] - round(chng, digits = 1)  
    } else if (Prchange > activeRow[2] && Prchange <= activeRow[3]) {
      forecastedYear[i+1,1] <- forecastedYear[i,1]
    } else if (Prchange > activeRow[3] && Prchange <= activeRow[4]) {
      chng <- runif(1,as.numeric(activeRow[7]),as.numeric(activeRow[8]))
      forecastedYear[i+1,1] <- forecastedYear[i,1] + round(chng, digits = 1)  
    } else {
      stop("Shit is fucked")
    }
  }
  
  ## then spring (an additional 91 days)
  for (i in 93:182) {
    if (i == 93 && forecastedYear[i,1] < 5.8) { forecastedYear[i,1] <- 5.8 }
    currentH <- forecastedYear[i,1]
    currentH <- round(currentH,digits = 1)
    activeRow <-  filter(altmarkSprDF, sprHeight == currentH) 
    Prchange <- runif(1,0,1)
    if (Prchange <= activeRow[2]) {
      chng <- runif(1,as.numeric(activeRow[5]),as.numeric(activeRow[6]))
      forecastedYear[i+1,1] <- forecastedYear[i,1] - round(chng, digits = 1) 
    } else if (Prchange > activeRow[2] && Prchange <= activeRow[3]) {
      forecastedYear[i+1,1] <- forecastedYear[i,1]
    } else if (Prchange > activeRow[3] && Prchange <= activeRow[4]) {
      chng <- runif(1,as.numeric(activeRow[7]),as.numeric(activeRow[8]))
      forecastedYear[i+1,1] <- forecastedYear[i,1] + round(chng, digits = 1) 
    } else {
      stop("Shit is fucked")
    }
  }  
  
 ## extract the 10th highest sillHeight for annualSeries
  forecastedYear <- forecastedYear*0.96
   annualSeries[ll,1] <- Rfast::nth(forecastedYear, 10, descending = T)
   annualSeries[ll,2] <- sum(forecastedYear[1:92,1])
  annualSeries[ll,3] <- sum(forecastedYear[93:182,1])
   
   forecastedYear[] <- NA
   print(ll)
}

write.csv(annualSeries,'C:/Workspace/MetapopHydrologyForecast.csv')




## This is going to be big
make sure to save : )

write.csv(annualSeries,'C:/Workspace/tsCenturies.csv')
## I have played with the data in excel ie multiply by 0.96 and cleaned a bit








