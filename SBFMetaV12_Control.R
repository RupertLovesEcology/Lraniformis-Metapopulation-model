##############################################################################################################
## The Metapopulation Model setup section
## This model will progress multiple discrete populations of different sizes with the capacity to move between sites each year
## We include 4 possible populations sizes (small, medium, large and extra-large)
## we include 3 levels of density feedback
## The wet dry sequence is forecasted using a modified Markov-chain for each of winter and spring to create annual sequences
## unsure of movement

Consider adding Wetness changes to movement 

Review  Survival during movement. Currently set to 0.5 survival

Model currently runs for 15 years of forecasted hydrology as a burn in
This is followed by 23 years of the Millenium drought hydrology
This decimates the population and creates 10000 starting scenarios for our subsequent treatments (found in PostDrought.csv  dim = 23*6,10000)
This version of the code takes those starting populations and runs them through the control scenario of 30 years of forecasted hydrology with no eWatering
final outputs are exported into finalControl.csv which is dim = 23*3,10000. The 23 * 3 is 1 column of each wetlands totalpopulation and 
then the adultpopulation of each wetland, then the total number of ewaters to eah wetland (will be 0 for the control treatment) 
Subsequent treatments will add  the number of eWater events for each wetland/treatment

WILL HAVE TO FIX WETtRACKER AND ADD A VARIABLE FOR TRACKING THE INTERVENTIONS FOR FUTURE TREATMENTS

Have commented out the multicore code for now
\


V10 i HAVE CHANGED THE OUPUT NUMBER AND  FORMATS BUT HAVE NOT REKNITTED IN THE OUTPUTS TO THE ARRAYS
I MAY ALSO NEED HOLDERS AND TO RESET STRATEGICALLY

V11 I have added patch quality to remove permanent wetlands from the privelidged position that drives the stability

##############################################################################################################

##################################################
## All of the requisite background code at the top
#################################################

## Rupert Mathwin & CJA Bradshaw
## Jan 2021


## Remove everything
rm(list = ls())



.libPaths("C:/workspace/math0286/R/win-library/4")
source("C:/workspace/math0286/R/win-library/3.6/matrixOperators.r")
setwd("C:/workspace/math0286/R")

## For the LPC1 computer
#source("C:/Workspace/RLibraries/4/matrixOperators.r")
#setwd("C:/Workspace/")
#.libPaths("C:/Workspace/RLibraries/4")

## libraries
library(DescTools)
library(pracma)
library(crayon)
library(DataCombine)
library(plyr)
library(adehabitatLT)
library(beepr)
library(readr)
library(parallel)
library(doSNOW)
library(iterators)
library(snow)
library(reshape2)



###############################################################################################################
## SET THE MODEL PARAMETERS

## movementTypes are A=DS along river, B=DS cross river, C=Normal overland, D=US along river,E=US cross river,F=50 km +
movementType <- read.csv("C:/Workspace/movementType.csv", header = FALSE, sep = ",", dec = ".")
wetlandMetadata <- read.csv("C:/Workspace/wetlandMetadataV2.csv", header = TRUE, sep = ",", dec = ".")
PostDrought <- read.csv("C:/Workspace/PostDrought.csv", header = TRUE, sep = ",", dec = ".")
## Sill heights are calculated for 600,000 years column 2 is winter wetness (MEAN=696) and column 3 is spring wetness (MEAN=760), column 4 is wint+spr wetness (MEAN=1447)
sillForecast <- read.csv("C:/Workspace/OrderedHydroForecast.csv", header = TRUE, sep = ",", dec = ".")

annualWetness <- read.csv("C:/Workspace/OrderedAnnualWetness.csv", header = TRUE, sep = ",", dec = ".")
# correct forecasted to match observed
annualWetness <- annualWetness * 0.98065
stayWet <- read.csv("C:/Workspace/stayWet.csv", header = T, sep = ",", dec = ".")

antecedentWet <- logical(length=23)


iter <- 100
generations <- 60
wetlandNum <- nrow(wetlandMetadata)

###  STORAGE FOR ERRORCHECKING (population size - S,M,L,XL)
runCounter <- array(data = 0, dim = wetlandNum)

## finalControl is the output.csv for this run
## 23*3 rows  23*1 is totalPop, 23*2 is adultPop, 23*3 is eWaterInterventions
#finalControl <- array(data = 0, dim = 23*3*10000)
#dim(finalControl) <- c(23*3,10000)
#class(finalControl) <- "numeric"

## Fudge factors for movement types A-G *note F is a reserved term for FALSE*
## movementTypes are A=DS along river, B=DS cross river, C=Normal overland, D=US along river,E=US cross river,G=50 km +
A <- 0.5
B <- 0.75
C <- 1
D <- 1.25
E <- 1.5
G <- 1

## RV of 17 spp. of anuran papers mean dispersal rates 16% +/- 14%
## Base probability of movement (will be adjusted)
wetMove <- 0.16

## mortality probability during movement
moveSurv <- 0.5

###############################################################################################################



#############################################################################################################
## Custom Functions Live Below
#############################################################################################################
# beta distribution shape parameter estimator function
##Generates an alpha and a beta value to inform the beta distribution 
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta)) }

##  Not used in this Southern Bell Frog Model
AICc.glm <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]}
  return(AICc.vec)}

## Not used in this Southern Bell Frog Model
k.glm <- function(x) {
  if (length(x$df.residual) == 0) k <- sum(x$dims$ncol) else k <- (length(x$coeff)+1)}

## Not used in this Southern Bell Frog Model
delta.IC <- function(x) x - min(x) ## where x is a vector of an IC

## Not used in this Southern Bell Frog Model
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC

## Not used in this Southern Bell Frog Model
chdev.glm <- function(x) ((( as.numeric(x[12]) - as.numeric(x[10]))/ as.numeric(x[12]))*100) ## % change in deviance, where x is glm object

# recursive function used when assigning to the immigrants holding array
SelectAge <- function(gens) {
  a <- round(runif(1,1,5))
  if (gens[a] > 0) {
    return(a)
  } else {
    SelectAge(emi.gen)
  }
}


## function to work out where each of the emigrants ends up
## use # destination <- determine.destination(??the place we leave from,emigrants,adj.wetlandMovement)
determine.destination <- function(depart,emigrnts,adjMove) { 
  ## Create a list of Levy flight distances for the moving frogs
  moveDist <- sample(dist.vec, emigrnts, replace=T, prob=pr.pred)
  # print(moveDist)
  ## WILL NEED TO ADD THE WETNESS/DISTANCE MODIFICATION BELOW
  # apologies to anyone reading for this next code snippet. a is the distance to the next DS site, b the distance to the next US site
  a <- 0
  b <- 0
  
  if (is.na(adjMove[depart,1])) {
    a <-  adjMove[depart,2] 
  } else if (is.na(adjMove[depart,2])) {
    a <- adjMove[depart,1] + abs(adjMove[depart,3] - adjMove[depart,1]) 
  } else {
    a <- adjMove[depart,1] + abs(adjMove[depart,2] - adjMove[depart,1])   
  }
  
  if (is.na(adjMove[depart,22])) {
    b <-  adjMove[depart,23] + abs(adjMove[depart,21] - adjMove[depart,23])   
  } else if (is.na(adjMove[depart,23])) {
    b <-  adjMove[depart,22] 
  } else {
    b <-  adjMove[depart,23] + abs(adjMove[depart,23] - adjMove[depart,22])  
  }
  
  ## c is the placeholder for losses upstream (mortality not applied), d is the placeholder for losses downstream 
  c <- 0
  d <- 0
  
  for (ll in 1:length(moveDist)) {
    #track if I need to 'next' (when the frog has left the reach)
    tracks <- 0
    current <- moveDist[ll]
    # cat("moveDist is ", moveDist, "\n")
    
    ## as I pass each distance I will check if it is larger than the distance to the downstream of upstream threshold
    # If yes then I will 50/50 check of the frog moves out of the reach (to be replaced at a different step) 
    if (current > a && runif(1,0,1) >= 0.5) {
      d <- d + 1
      tracks <- 1
      #  cat("frog ", ll, " was lost from wetland ", wetlands, " downstream \n")
    } 
    if (tracks == 0 && current > b && runif(1,0,1) >= 0.5) {
      c <- c + 1
      tracks <- 1
      #  cat("frog ", ll, " was lost from wetland ", wetlands, " upstream \n")
    } 
    if (tracks == 1) {
      moveDist[ll] <- NA
      next
    }
    
    tracks <- 0
    smallest <- 50000
    equiDist <- 0
    #  cat("Starting frog number ", ll, "\n")
    #  cat("moveDist is ", moveDist, "\n")
    for (mm in 1:ncol(adjMove)) {
      #  cat("Checking column ", mm, "\n")
      ddd <- abs(current - adjMove[depart,mm])
      #  cat("ddd is ", ddd, "\n")
      if (ddd < smallest && !is.na(ddd)) {
        smallest <- ddd
        moveDist[ll] <- mm
        #   cat("replacing ddd moveDist position ", ll, " with the chosen column, which is ", mm, "\n")
      } }
    # check for multiple sites with the same minimum distance
    for (nn in 1:ncol(adjMove)) {
      #  cat("checking for equal values \n")
      if (!is.na(adjMove[depart,nn]) && adjMove[depart,nn] == smallest) { 
        equiDist <- cbind(equiDist,nn)
      } 
      ## check all for equal values then randomise one of them assign the final destination to moveDist
      #    cat("equiDist is ", equiDist, "\n")
      if (length(equiDist) > 2) {
        chosen <- floor(runif(1, min=2, max=ncol(equiDist))) 
        #   cat("I have selected column", chosen, " this corresponds to a value of ", equiDist[chosen],"\n")
        moveDist[ll] <- equiDist[chosen] 
      }  
    }
    
  } 
  # cat("c is ", c, "\n")
  # cat("d is ", d, "\n")
  # cat("movedist is \n")
  # print(moveDist)
  
  # quick reminder c is US loss, d is DS loss
  moveDist <- as.numeric(moveDist[!is.na(moveDist)])
  moveDist <- append(moveDist,c)
  moveDist <- append(moveDist,d)
  # cat("movedist is \n")
  # print(moveDist)
  return(moveDist)
}


## Assign immigrants from emi.gen into immigrants, calculating movement survival on the way past
#use "immigrants <- add.immigrants(emi.gen,moveSurv,destinations,immigrants) "
add.immigrants <- function(emigrntGen,survival,desti,immiArray,storeMove,storeRecol,it) { 
  #  cat("starting add.Immigrants function \n")
  for (oo in 1:length(emigrntGen)) {
    #  cat("oo is ", oo, "\n")
    while (emigrntGen[oo] > 0) {
      #  cat("taking one \n")
      emigrntGen[oo] <- emigrntGen[oo] - 1
      # check if it survived
      if (runif(1,0,1) <= survival) {
        #  cat("survived")
        immiArray[desti[1],oo] <- immiArray[desti[1],oo] + 1
        desti <- desti[-1]
      }
    }
  }
  #  cat("finishing add.Immigrants function \n")
  return(immiArray)
}


## determine age of immigrants DetermineAges(emiLoss)
#code snippet to assign ageClass according to expected age structure (noting range 1 is 1-2yrs, 2 is 2-3 yrs)
DetermineAges <- function(immiNums) {
  #  cat("starting DetermineAges \n")
  fromUS <- NA
  fromDS <- NA
  counter <- 1 
  for (sec in 1:2) {
    while (counter <= immiNums[1,sec]) {
      ## THIS STEP ADDS THE MORTALITY DURING MOVEMENT NEED TO CHANGE ONCE WE FIGURE OUT WHAT THAT LOOKS LIKE
      if (runif(1,0,1) < 0.5) {
        counter <- counter + 1
        next 
      }
      ageImi <- runif(1,0,1)
      if (ageImi <= 0.66) {
        ifelse(sec == 1,fromDS <- cbind(fromDS,1),fromUS <- cbind(fromUS,1))
      } else if (ageImi > 0.66 && ageImi <= 0.93) {
        ifelse(sec == 1,fromDS <- cbind(fromDS,2),fromUS <- cbind(fromUS,2))
      } else if (ageImi > 0.93 && ageImi <= 0.987) {
        ifelse(sec == 1,fromDS <- cbind(fromDS,3),fromUS <- cbind(fromUS,3))
      } else {
        ifelse(sec == 1,fromDS <- cbind(fromDS,4),fromUS <- cbind(fromUS,4))
      }
      counter <- counter + 1
    }
    counter <- 1 
  }
  fromDS <- fromDS[-1]
  fromUS <- fromUS[-1]
  #  cat("fromDS is ",fromDS,"\n")
  #  cat("fromUS is ",fromUS,"\n")
  
  if (length(fromUS) > length(fromDS)) {
    immiAges <- array(data=NA, dim = 2 * length(fromUS))
    dim(immiAges) <- c(2,length(fromUS))
    rownames(immiAges) <- c("from DS", "from US")
  } else {
    immiAges <- array(data=NA, dim = 2 * length(fromDS))
    dim(immiAges) <- c(2,length(fromDS))
    rownames(immiAges) <- c("from DS", "from US")
  }
  i <- 1
  while (i <= length(fromDS)) {
    immiAges[1,i] <- fromDS[i]
    i <- i + 1
  } 
  i <- 1
  while (i <= length(fromUS)) {
    immiAges[2,i] <- fromUS[i]
    i <- i + 1
  } 
  # cat("finishing DetermineAges \n") 
  return(immiAges)
}




## set time limit for projection in 1-yr increments
yr.now <- 2020
#************************
#yr.end <- 2020 + round((40*gen.l), 0) # set projection end date

yr.end <- 2020 + 86
#************************
##Note I have constrained t to 86 years
t <- (yr.end - yr.now)
yrs <- seq(yr.now,yr.end,1)
longev <- 5
age.vec <- seq(0,longev,1)
lage <- length(age.vec)
sex.ratio <- 0.5
stages <- lage
## set population storage matrices n.mat
n.mat <- array(data = 0, dim = (stages * (t + 1) * wetlandNum))
dim(n.mat) <- c(stages,(t + 1),wetlandNum)
popmat <- matrix(0,nrow=stages,ncol=stages)
colnames(popmat) <- age.vec[1:stages]
#rownames(popmat) <- age.vec[1:stages]
rownames(popmat) <- c("Fecundity","Survival to 1 year","Survival to 2 year","Survival to 3 year","Survival to 4 year","Survival to 5 year")


## fertility data 
# clutch.size.la <- c(4124,6178) # L. aurea
clutch.size.lr <- c(1885,3893) # L. raniformis
clutch.size.lr.sd.prop <- 502/mean(c(1885,3893))




prop.breeding <- c(0,rep(1,5))
fert.mn <- mean(clutch.size.lr)*prop.breeding
fert.mn

#duration data (eggs and tadpoles) following line is for the LHC
hatch.dur <- c(2,4)
hatch.dur.sd.prop <- 0.5/3
tadpole.dur <- c(70,80) # duration (days) L. raniformis 23 deg
tadpole.dur.sd.prop <- 2.5/75

## survival data
hatch.pr <- c(0.933, 1) # 2-4 days, L. raniformis
tadpole.s.cs <- 0.10 # C. signifera
tadpole.s.ra <- 0.05 # R. aurora
ad.s.daily <- c(0.975, 0.996) # daily survival L. raniformis

##the below uses the more complete survival TO METAMORPHOSIS figures from Bull (C.signifera)
tadpole.mn.1 <- mean(c(.15,.26))
tadpole.sd.1 <- ((.26-tadpole.mn.1)+(tadpole.mn.1-.15))/2/1.96
tadpole.mn.2 <- mean(c(.07,.56))
tadpole.sd.2 <- ((.56-tadpole.mn.1)+(tadpole.mn.1-.07))/2/1.96
tadpole.mn <- mean(c(tadpole.mn.1, tadpole.mn.2))
tadpole.sd <- sqrt(tadpole.sd.1^2 + tadpole.sd.2^2)
tp.s.alpha <- estBetaParams(tadpole.mn, tadpole.sd/10)$alpha
tp.s.beta <- estBetaParams(tadpole.mn, tadpole.sd/10)$beta

tomet.dur.iter <- round(sum(c(runif(1, hatch.dur[1], hatch.dur[2]), runif(1, tadpole.dur[1], tadpole.dur[2]))), 0)
toad.dur.iter <- 365 - tomet.dur.iter
#tomet.s.iter <- (runif(1, min=hatch.pr[1], max=hatch.pr[2])) * (runif(1, min=tadpole.s.ra, max=tadpole.s.cs)) 

tomet.dur.mn <- round(sum(c(mean(hatch.dur), mean(tadpole.dur))), 0)
toad.dur.mn <- 365 - tomet.dur.mn

#tomet.s.mn <- mean(hatch.pr) * mean(c(tadpole.s.ra, tadpole.s.cs))
#juv.s.daily.vec <- runif(toad.dur.iter, min=ad.s.daily[1], max=ad.s.daily[2])
#juv.s.daily.mn <- rep(mean(ad.s.daily), toad.dur.mn)
#ad.s.iter <- prod(runif(365, min=ad.s.daily[1], max=ad.s.daily[2]))
#ad.s.mn <- (mean(ad.s.daily))^365

#NEW ADULT ANNUAL SURVIVAL VALUES FROM PICKETT L. AUREA
ad.s.yr.mn <- 0.2172 # Litoria aurea (Pickett et al.)
ad.s.yr.sd <- 0.087 # Litoria aurea (Pickett et al.)

##FROM PICKETT
ad.s.yr.alpha <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$alpha
ad.s.yr.beta <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$beta

##Calculate our survivals

#likelyhood of egg hatching and tadpole surviving to metamorphosis
tomet.s.iter <- (rbeta(1, tp.s.alpha, tp.s.beta)) * (runif(1, min=0.933, max=1))

#sample a daily probability of survival
toad.daily.s.iter <- nthroot(rbeta(1,ad.s.yr.alpha, ad.s.yr.beta) , 365)

toad.s.season.iter <- toad.daily.s.iter ^ toad.dur.mn
## change the above line to the below line just need the syntax to resample
## toad.s.season.iter <- prod(runif(toad.dur.mn, toad.daily.s.iter))
toad.s.iter <- tomet.s.iter * toad.s.season.iter                     

# TO RESAMPLE annual survival of an adult
#  ad.s.iter <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)

#Create the survival vector, to adult survival then 4 adult survivals
ad.s.vec.iter <- rep(NA,5)
for (s in 1:5) {
  ad.s.vec.iter[s] <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)  }
surv.iter <- c(toad.s.iter, ad.s.vec.iter)



# We hereafter consider 4 population sizes, each represented with 4 density feedback values (one for each of egg, tad, juve, adult)
# Small, medium, large, very large.
# Assign the number of spawning masses that each site holds and convert it to the number of female egges (uses mean clutch size)
K.rel.egg.dem.vec <- c(50, 150, 500, 1000)
K.rel.egg.dem.vec <- K.rel.egg.dem.vec * 1444

# tadpoles surviving to 1 year old will inhibit themselves NOTE THIS EQUATION ONLY CONSIDERS FEMALES SO THE ACTUAL NUMBER WILL BE DOUBLE  
K.rel.tad.dem.vec <- c(100, 300, 800, 1000)

# number of first years present will influence survival from 1 to 2 years old but more strongly than for older age brackets         
K.rel.juv.dem.vec <- c(200, 600, 1800, 3000)

# number of first years present will influence adult survival NOTE THIS EQUATION ONLY CONSIDERS FEMALES SO ACTUAL NUMBER WILL BE DOUBLE 
#These numbers are used for the probability of emigration (ie above K high movement, below K lower movement) 
K.rel.adult.dem.vec <- c(5000,1000,6000,11550)

initSpawn <- c(10,30,85,150)

#set founding fem pop sizes 
init.fem.pop <- c(20,65,150,500) 

init.vec <- rep(NA,24)
dim(init.vec) <- c(6,4)

##Populate the Matix (popmat) and create a failure matrix(popmat.fail) for years with no breeding
## removed  diag(popmat[2:(stages), ]) <- surv.mn[-stages]
diag(popmat[2:(stages), ]) <- surv.iter[-stages]
popmat[stages,stages] <- 0 # surv.mn[stages] 
popmat[1,] <- fert.mn * sex.ratio
popmat.orig <- popmat ## save original matrix as popmat.orig
popmat <- popmat.orig



#Create init.vec note 2nd dimension is 1 = S, 2 = M, 3 = L, 4 = XL
## This will start every wetland with the same size-appropriate starting populations.
## I can recalculate each init.vec if necessary
totalSurv <- sum(popmat[2:6,1:6])
for (size in 1:4) {
  for (x in 1:6) { init.vec[x,size] <- 0 }
  init.vec[1,size] <- round(initSpawn[size] * (runif(1, min=clutch.size.lr[1], max=clutch.size.lr[2])))
  for (dd in 2:6) { init.vec[dd,size] <- round((sum(popmat[2:6,(dd - 1)])/totalSurv) * (init.fem.pop[size]))    } }

## Assign the starting population for each of the 23 wetlands as determine by their size category in n.mat
#for (size in 1:wetlandNum) { 
#  n.mat[,1,size] <- init.vec[,wetlandMetadata[size,1]] 
#}


popmat.current <- popmat
popmat.fail <- popmat
popmat.fail[1,] <- 0 

#create the egg density correcting function and variables
#   line equation y = 1.01 - ((x ^ (9/x)) * (x ^ 2)))
surv.mult.egg.up <- 1.0
surv.mult.egg.upmid <- 1.0
surv.mult.egg.mid <- 0.94
surv.mult.egg.lo <- 0.8
surv.mult.egg.lo.lo <- 0.2
surv.mult.egg.lo.lo.lo <- 0.05

K.egg.up <- 1
K.egg.upmid <- 0.98
K.egg.mid <- 0.90
K.egg.lo <- 0.7
K.egg.lo.lo <- 0.3
K.egg.lo.lo.lo <- 0.01

K.egg.vec <- c(K.egg.up,K.egg.upmid, K.egg.mid,K.egg.lo, K.egg.lo.lo, K.egg.lo.lo.lo)
surv.mult.egg.vec <- rev(c(surv.mult.egg.up, surv.mult.egg.upmid, surv.mult.egg.mid, surv.mult.egg.lo, surv.mult.egg.lo.lo, surv.mult.egg.lo.lo.lo))
plot(K.egg.vec, surv.mult.egg.vec, pch=19)
##plot(K.egg.vec, surv.mult.egg.vec, type="n")
DD.dat <- data.frame(K.egg.vec, surv.mult.egg.vec)
param.init <- c(1.01, 9, 2)
fit.expd.egg <- nls(surv.mult.egg.vec ~ (a - ((K.egg.vec ^ (b/K.egg.vec)) * (K.egg.vec ^ d))),
                    algorithm = "port",
                    start = c(a = param.init[1], b = param.init[2], d = param.init[3]),
                    trace = TRUE,      
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
K.pred.egg.vec <- seq(0.01, 1, 0.01)
pred.surv.egg.mult <- (as.numeric(coef(fit.expd.egg)[1]) - (K.pred.egg.vec ^ (as.numeric(coef(fit.expd.egg)[2])/K.pred.egg.vec) * (K.pred.egg.vec ^ (coef(fit.expd.egg)[3]))))
lines(K.pred.egg.vec, pred.surv.egg.mult, lty=2, lwd=1, col="red")

eggfunc <- function(x) (1.01 - ((x ^ (9/x)) * (x ^ 2)))
intArea <- integrate(eggfunc,lower=0,upper=0.99)
area1 <- (as.numeric(intArea[1])/0.99)
s.mult.egg.iter <- (as.numeric(coef(fit.expd.egg)[1]) - (0.999 ^ (as.numeric(coef(fit.expd.egg)[2])/0.999) * (0.999 ^ (coef(fit.expd.egg)[3]))))

# plot(K.egg.vec, surv.mult.egg.vec, type="n", xlab="Carrying capacity of eggs(K)", ylab="Reduction in egg laying success")
# lines(K.pred.egg.vec, pred.surv.egg.mult, lty=1, lwd=1, col="black")

# plot(K.juv.vec, surv.mult.juv.vec, type="n", xlab="Proportion of the site carrying capacity (K)", ylab="Reduction in survival rate")
# lines(K.pred.juv.vec, pred.surv.juv.mult, lty=1, lwd=1, col="black")

## invoke a density-feedback function on tadpole survival to year 1
# density feedback survival multiplier for tadpoles hinges on the density of other tadpoles in the pond
surv.mult.up <- 1.0
surv.mult.upmid <- 0.58
surv.mult.mid <- 0.19
surv.mult.lo <- 0.10

K.up <- 1
K.upmid <- 0.83
K.mid <- 0.45
K.lo <- 0.3

K.tad.vec <- c(K.up,K.upmid, K.mid,K.lo)
surv.mult.tad.vec <- rev(c(surv.mult.up, surv.mult.upmid, surv.mult.mid, surv.mult.lo))
# plot(K.tad.vec, surv.mult.tad.vec, pch=19)

# Bleasdale
# y = (a + bx)^(-1/c)
DD.dat <- data.frame(K.tad.vec, surv.mult.tad.vec)
param.init <- c(-2.41e-01, 1.54, 1.17)
fit.expd.tad <- nls(surv.mult.tad.vec ~ (a + (b*K.tad.vec))^(-1/c), 
                    data = DD.dat,
                    algorithm = "port",
                    start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                    trace = TRUE,      
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
# plot(K.tad.vec, surv.mult.tad.vec, pch=19, xlab="K", ylab="reduction Tadpole survival to 1 yr")
K.pred.tad.vec <- seq(K.lo,1,0.01)
pred.surv.tad.mult <- (as.numeric(coef(fit.expd.tad)[1]) + (K.pred.tad.vec * as.numeric(coef(fit.expd.tad)[2])))^(-1/as.numeric(coef(fit.expd.tad)[3]))
# lines(K.pred.tad.vec, pred.surv.tad.mult, lty=2, lwd=1, col="red")

## invoke a density-feedback function on Juvenile survival from year 1 to year 2
# density feedback survival multiplier for juveniles hinges on the density of themselves (but more strongly than adults)
surv.mult.up <- 1.0
surv.mult.upmid <- 0.58
surv.mult.mid <- 0.19
surv.mult.lo <- 0.10

K.up <- 1
K.upmid <- 0.83
K.mid <- 0.45
K.lo <- 0.3

K.juv.vec <- c(K.up,K.upmid, K.mid,K.lo)
surv.mult.juv.vec <- rev(c(surv.mult.up, surv.mult.upmid, surv.mult.mid, surv.mult.lo))
# plot(K.juv.vec, surv.mult.juv.vec, pch=19)

# Bleasdale
# y = (a + bx)^(-1/c)
DD.dat <- data.frame(K.juv.vec, surv.mult.juv.vec)
param.init <- c(-2.41e-01, 1.54, 1.17)
fit.expd.juv <- nls(surv.mult.juv.vec ~ (a + (b*K.juv.vec))^(-1/c), 
                    data = DD.dat,
                    algorithm = "port",
                    start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                    trace = TRUE,      
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
#  plot(K.juv.vec, surv.mult.juv.vec, pch=19, xlab="Proportion of the site carrying capacity (K)", ylab="Reduction in survival rate")
#  plot(K.juv.vec, surv.mult.juv.vec, type="n", xlab="Proportion of the site carrying capacity (K)", ylab="Reduction in survival rate")
K.pred.juv.vec <- seq(K.lo,1,0.01)
pred.surv.juv.mult <- (as.numeric(coef(fit.expd.juv)[1]) + (K.pred.juv.vec * as.numeric(coef(fit.expd.juv)[2])))^(-1/as.numeric(coef(fit.expd.juv)[3]))
#  lines(K.pred.juv.vec, pred.surv.juv.mult, lty=2, lwd=1, col="red")
#  plot(K.juv.vec, surv.mult.juv.vec, type="n", xlab="Proportion of the site carrying capacity (K)", ylab="Reduction in survival rate")
#  lines(K.pred.juv.vec, pred.surv.juv.mult, lty=1, lwd=1, col="black")

fit.expd.adult <- fit.expd.tad


## Carrying capacity is determined by Pr(move) = n(t)/K^2 * 30 *note this is an individual probability of movement, not a proportion of movement
# where K.rel.adult.dem.vec[size(1-4)] is K
# Calculated within the code loop


## Levy flight Setup code
## use 'dist.rand1 <- sample(dist.vec, 100, replace=T, prob=pr.pred)' to produce 100 values
## to create dist.rand1

min.dist <- 50
max.dist <- 50000
lmax.dist <- log(max.dist)
max.prob.mov <- 0.3
nmoves <- 100000
u <- simm.levy(1:nmoves, l0=min.dist, mu = 2.3)
hist.out <- hist(log(u[[1]]$dist),main="", xlab="log distance (m)", ylab="frequency", br=50)
range(u[[1]]$dist, na.rm=T)

lcounts <- log(hist.out$counts)
lprob.scaled <- log(max.prob.mov * hist.out$counts/max(hist.out$counts))
ldist <- hist.out$mids
ldt <- na.omit(data.frame(ldist,lcounts,lprob.scaled))
ldat.clean <- ldt[which(ldt$ldist <= lmax.dist), ]

if (length(which(is.infinite(ldat.clean$lcounts)==T)) > 0) {
  ldat.rem <- which(is.infinite(ldat.clean$lcounts)==T)
  ldat.clean <- ldat.clean[-ldat.rem,]
}

plot(ldat.clean$ldist, ldat.clean$lprob.scaled, pch=19, ylab="log scaled prob", xlab="log dist (m)")
fit.levy <- lm(lprob.scaled ~ ldist, data=ldat.clean)
abline(fit.levy, lty=2, col="red")
summary(fit.levy)
exp(range(ldat.clean$ldist, na.rm=T))

levy.int <- as.numeric(coef(fit.levy)[1])
levy.sl <-  as.numeric(coef(fit.levy)[2])

dist.vec <- seq(min.dist, max.dist, 1)
pr.pred <- exp(levy.int + (log(dist.vec))*levy.sl)
plot(dist.vec, pr.pred, type="l")

dist.rand <- sample(dist.vec, 1, replace=T, prob=pr.pred)
dist.rand
#hist(dist.rand)



######
## Populate underlying movement data (site-to-site distance which incorporates site-to-site movement type)
######
wetlandMovement <- read.csv("C:/Workspace/wetlandMovement.csv", header = FALSE, sep = ",", dec = ".")
## And a temporary holder for the seasonally adjusted matrix 
bass <- seq(1,nrow(wetlandMovement),1)
bass1 <- bass2 <- bass

for (i in 1:wetlandNum) {
  bass1[i] <- paste("From",bass[i])
}
rownames(wetlandMovement) <- bass1

for (i in 1:wetlandNum) {
  bass2[i] <- paste("To",bass[i])
}
colnames(wetlandMovement) <- bass2
adj.wetlandMovement <- wetlandMovement
adj.wetlandMovement[] <- NA

### Create an array of distance between wetlands where the distance is adjusted to reflect the type of journey ie downstream easier/shorter, upstream harder/further.

for (i in 1:nrow(wetlandMovement)) {
  for (j in 1:ncol(wetlandMovement)) {
    if (i == j) {
      adj.wetlandMovement[i,j] <- NA
    } else {
      adj.wetlandMovement[i,j] <- wetlandMovement[i,j]*get(movementType[i,j])
      if (adj.wetlandMovement[i,j] > 50) {  adj.wetlandMovement[i,j] <- 50 }
    }
  }
}

## and convert to metres
adj.wetlandMovement <- adj.wetlandMovement * 1000

## Create a holder for losses out the top and bottom of the reach during emmigration
emiLoss <- array(data = 0, dim = 2)
dim(emiLoss) <- c(1,2)
colnames(emiLoss) <- c("DS Loss","US Loss")

exoImmi <- array(data = NA, dim = 2)
dim(exoImmi) <- c(1,2)
colnames(exoImmi) <- c("Wetland","AgeClass")

## Create an Array to hold all of the successful movements
## for the following arrays we store the cumulative movements/recolonisations/ewater/adult population
## at the end of each decade - 10,20, . . . 60
successfulMovement <- array(data = 0, dim = (23*23*6))
dim(successfulMovement) <- c(23,23,6)
successfulMovementHold <- array(data = 0, dim = 23*23)
dim(successfulMovementHold) <- c(23,23)

## Repeat this array for all of the successful movements that recolonised a locally extinct wetland
successfulRecolonise <- array(data = 0, dim = (23*23*6))
dim(successfulRecolonise) <- c(23,23,6)
successfulRecoloniseHold <- array(data = 0, dim = 23*23)
dim(successfulRecoloniseHold) <- c(23,23)

## Repeat this array for all of the individual patch extinction events
patchExtAll <- array(data = 0, dim = (23*7))
dim(patchExtAll) <- c(23,7)
patchExt <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)


## and one for eWater counting
eWater <- array(data = 0, dim = (23*6))
dim(eWater) <- c(23,6)
eWaterHold <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

## AdultPop and AllPop will have an additional row 
## Row24 is the cumulative total of how many pops were alive at each 10yr step
## Create an array to hold all of the final populations at end of each iter/10
adultPop <- array(data = 0, dim = 24 * 7)
dim(adultPop) <- c(24,7)

allPop <- array(data = 0, dim = 24 * 7)
dim(allPop) <- c(24,7)

## These next arrays will manage the eWater treatments
## Create dryTracker to hold the wetness/dryness for each wetland. We use this to decide on interventions 
#dryTracker <- 1:23
#dryTracker[] <- 0

## SDArray for calculating StDev/CIs
# 3rd dimensions - 1 = adult pop, 2 = eWater, 3 patch extinctions, 4 all movement, 5 all recolonisations,
SDArray <- array(data = 0, dim = 5 * 10000 * 7)
dim(SDArray) <- c(5,7,10000)

# SDWetlandMove will track the number of movements and recolonisations from each wetland 
# dims are 23wetlands,7 time steps,2 (1move,2recol),10000 iters
SDWetlandMove <- array(data = 0, dim = 23*7*2 * 10000)
dim(SDWetlandMove) <- c(23,7,2,10000)

wetCycle <- array(data=0,dim=23*5)
dim(wetCycle) <- c(5,23)

#eWaterTracker <- 1:5
#eWaterTracker[] <- 0

eWaterSites <- c(0)

#Spread eWatering sites
#  eWaterSites <- c(4,8,13,21)

#Clumped eWatering sites
# eWaterSites <- c(11,13,14,16)
 
 banrockFish <- c(F,T,T)
 
 
 ## Array to hold for probability of occupancy ?heatmap
 occupancyHeatC <- array(data=0,dim=23*7*10000)
 dim(occupancyHeatC) <- c(23,7,10000)
 

 ## TO USE THE above 'wetlands %in% eWaterSites' returns T/F
## match(wetlands,eWaterSites) returns the index location ie to feed into eWaterTracker


## The probability of mortality whilst moving between sites
## NOTE  THIS MORTALITY IS APPLIED ABOVE AND BEYOND BACKGROUND MORTALITY 
## ie NORMAL MORTALITY HAS ALREADY BEEN APPLIED AND THIS MOVEMENT MORTALITY IS APPLIED ON TOP OF THIS WITHIN THE SAME YEAR 
## CURRENTLY USING A SINGLE VARIABLE moveSurv
#if I want to revert back to an array it is below
#movementSurv <- array(data = 0.55,dim = wetlandNum*wetlandNum)
#dim(movementSurv) <- c(wetlandNum,wetlandNum)
#rownames(movementSurv) <- bass1
#colnames(movementSurv) <- bass2









###############################################################################################################
## And NOW for the actual model!
###############################################################################################################



## Version 4.0
##############################################################
## 23 populations of four possible sizes (S,M,L,XL)
## 60 generations 50 + 10 for burn in
### UPDATE THIS SECTION ON COMPLETION
## sill threshold for the year is sourced from sillHeight (600,000 calculated sill heights and sum of winter [col1] and sum of spring[col2] wetness)
## adj.wetlandMovement stores the modified distance(actual distance adjusted by movement type e.g. overland(harder),downstream(easier) variables A-G)
## function 'determine.Destination' randomises the length of each emigration, determines the potential destinations and then randomly selects a destination
## 0.55 Pr of surviving in transit to all populations (stored in moveSurv)
## wetland names/sizes in wetlandMetadata

##############################################################
## To Do - Assign the correct level of density dependant inhibition !!!!!!!!!!!!!
## To Do - determine how to run the hyperspace cube to run the sensitivity analysis 15 - 25 variables

############################################################



## set up arrays to track emigration and immigration
moveAgeNm <- c("1-2yr","2-3yr","3-4yr","4-5yr","5+yr") 

## emigrants is from wetland# to wetland # and third dimension is age class (only 1-5 year olds move)
emigrants <- array(data = 0, dim = wetlandNum * wetlandNum * 5)
dim(emigrants) <- c(wetlandNum, wetlandNum, 5) 

emi.gen <- array(data = 0, dim = 5)
dim(emi.gen) <- c(1*5)

immigrants <- array(data = 0, dim = wetlandNum * 5)
dim(immigrants) <- c(wetlandNum,5)
rownames(immigrants) <- bass2
colnames(immigrants) <- moveAgeNm

## Set up parallel processing (nproc is the number of processing cores to use)
#nproc <- 6
#cl.tmp = makeCluster(rep('localhost', nproc), type='SOCK')
#registerDoSNOW(cl.tmp)
#getDoParWorkers()

SDTracker <- 0

immigrantDestinations <- 0

## The Outermost Loop: iterate the process iter times (spoiler alert it's one time)
## note I currently have 150 unique centuries if it exceeds 150 we will have to start again
for (e in 1:iter) {
  cat("starting iteration ", e, "\n")  
  SDTracker <- SDTracker + 1
  
  sillCounter <- 0 
  
  
  # reset the population matrix except for the starting populations (n.mat[,1,])
  # n.mat[1:6,2:ncol(n.mat),1:wetlandNum] <- 0 
  
  
  
  
  #reset the dryTracker
 # dryTracker[] <- 0 
  successfulMovementHold[] <- 0
  successfulRecoloniseHold[] <- 0
  eWaterHold[] <- 0
  patchExt[] <- 0
  antecedentWet[] <- F
  wetCycle[] <- 0
  
  
  ## The Second Loop: run the current projection set up for the number of generations (years actually) + 10 for burn-in 
  ## NOTE I have generated centuries so this value can go up to 90 yrs as required, currently at 50
  for (i in 1:(generations)) {
    ## work out the wetness figures for this year the sillHeight and wetMod for this year
    # note I have applied the reduction of 10cm to allow for the threshold in this step
    
    sillCounter <- sillCounter + 1
    sillHeight <- sillForecast[e,sillCounter]
   
    
    #    for (cycle in 1:wetlandNum) {
    #      BottomVal <- cycle * 6
    #      #  cat(" Bottomval is ", BottomVal, " n.mat[1:6,i,wetlands] is  ", n.mat[1:6,i,wetlands], "\n")
    #      PostDrought[(BottomVal - 5):BottomVal,e] <- n.mat[1:6,i,cycle]
    
    
    
    
    
    
    
    # cat(" generations is ", i, "\n")
    if (i == 1) {
      for (st in 1:wetlandNum) {
        BottomVal <- st * 6
        n.mat[1:6,i,st] <-  PostDrought[(BottomVal - 5):BottomVal,e]
        class(n.mat) <- "numeric"
        # add to adult
        adultPop[st,1] <- adultPop[st,1] + sum(n.mat[2:6,i,st])
        if (sum(n.mat[2:6,i,st]) > 0) { adultPop[24,1] <- adultPop[24,1] + 1 }
        
        ## add to all
        allPop[st,1] <- allPop[st,1] + sum(n.mat[1:6,i,st])
        if (sum(n.mat[1:6,i,st]) > 0) { allPop[24,1] <- allPop[24,1] + 1 }
        
        
        
        #  cat("n.mat starting pop is /n", n.mat[,1,], "\n")
              }
    }
    # cat("n.mat starting pop is /n", n.mat[,1,], "\n")  
    
    #   cat(red("starting generation ", i, "the sillHeight this year is ", sillHeight,"\n"))
    #iterate the runCounter used for errorchecking
    if (i==1) { runCounter[wetlandNum] <- (runCounter[wetlandNum] + 1) 
    }
    
    # if there are no frogs alive break from this loop
    if (sum(n.mat[,i,]) == 0) {
      cat("i is ", i, "  sillHeight is ", sillHeight, "  \n")
      print(red("All populations are extinct, which is sad"))
    #  beep(9)
      break  
    }
    
    
    
    
    
    ## Cycle through each of the wetlands, The Wetlands Loop
    for (wetlands in 1:wetlandNum) {
      
      
      
      #   cat("I am starting wetland ", wetlands, " which is size ", wetlandMetadata[wetlands,1], "\n")
      
      # resample the duration of egg to tadpole
      # tomet.dur.iter <- round(sum(c(runif(1, hatch.dur[1], hatch.dur[2]), runif(1, tadpole.dur[1], tadpole.dur[2]))), 0)
      tomet.dur.iter <- round((rnorm(1, hatch.dur, hatch.dur.sd.prop*hatch.dur) + rnorm(1, tadpole.dur, tadpole.dur.sd.prop*tadpole.dur)), 0)
      toad.dur.iter <- 365 - tomet.dur.iter
      
      ##Calculate our survivals
      #likelyhood of egg hatching and tadpole surviving to metamorphosis
      tomet.s.iter <- (rbeta(1, tp.s.alpha, tp.s.beta)) * (runif(1, min=0.933, max=1))
      toad.daily.s.iter <- nthroot(rbeta(1,ad.s.yr.alpha, ad.s.yr.beta) , 365)
      toad.s.season.iter <- toad.daily.s.iter ^ toad.dur.iter
      toad.s.iter <- tomet.s.iter * toad.s.season.iter                     
      ad.s.iter <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)
      
      #Create the survival vector (popmat) for the year (density dependence not considered)
      ad.s.vec.iter <- rep(NA,5)
      for (s in 1:5) {
        ad.s.vec.iter[s] <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)  }
      surv.iter <- c(toad.s.iter, ad.s.vec.iter)
      
      min.clutch.size.lr <- as.numeric(quantile(rnorm(n=1000, mean=clutch.size.lr, sd=clutch.size.lr.sd.prop*clutch.size.lr), probs=0.025, na.rm=T))
      max.clutch.size.lr <- as.numeric(quantile(rnorm(n=1000, mean=clutch.size.lr, sd=clutch.size.lr.sd.prop*clutch.size.lr), probs=0.975, na.rm=T))
     
      
      #  fert.iter <- round((runif(stages, min=clutch.size.lr[1], max=clutch.size.lr[2])) * prop.breeding, 0)
      # may have to change maxAge to stages
      fert.iter <- round((runif(stages, min=min.clutch.size.lr, max=max.clutch.size.lr)) * prop.breeding, 0)
      popmat[1,] <- fert.iter * sex.ratio 
      diag(popmat[2:(stages), ]) <- surv.iter[-stages]
      
      #note these are placeholder values for popmat.current, they are recalculated below not necesary but retained for safety
      popmat.current <- popmat
      
      
      
      # Implement density dependence effects for each of four densities and feed into the popmat.all[,,]
      #note density effect on egg laying is applied after matrix multiplication 
      #  density feedback for tadpoles to first year is strong and driven by the number of tadpoles in the cohort
      s.mult.iter.tad <- 1
      s.mult.iter.juv <- 1
      s.mult.iter.ad <- 1
      #   cat("i is ", i, "wetlands is", wetlands,"\n")
      
      # instil density dependence for tadpoles growing into year 1 adults
      K.rel.tad <- (n.mat[1,i,wetlands]/K.rel.tad.dem.vec[wetlandMetadata[wetlands,1]])
      #   cat(" K.rel.tad is ", K.rel.tad, "\n")
      if (!is.nan(K.rel.tad)) {
        if (K.rel.tad > 2.1)  { K.rel.tad <- 2.1 }
        if (K.rel.tad <= 2.1) {
          s.mult.iter.tad <- (as.numeric(coef(fit.expd.tad)[1]) + (K.rel.tad * as.numeric(coef(fit.expd.tad)[2])))^(-1/as.numeric(coef(fit.expd.tad)[3]))
          popmat.current[2,1] <- popmat[2,1] * s.mult.iter.tad   } 
      } 
      
      # instil density dependence for juveniles (1 - 2 years) is driven by the  number of yr 1 present/competing per 
      K.rel.juv <- (n.mat[2,i,wetlands]/K.rel.juv.dem.vec[wetlandMetadata[wetlands,1]]) 
      if (!is.nan(K.rel.juv)) {
        
        if (K.rel.juv > 2.1)  { K.rel.juv <- 2.1 }
        
        if (K.rel.juv <= 2.1) {
          s.mult.iter.juv <- (as.numeric(coef(fit.expd.juv)[1]) + (K.rel.juv * as.numeric(coef(fit.expd.juv)[2])))^(-1/as.numeric(coef(fit.expd.juv)[3]))
          popmat.current[3,2] <- popmat[3,2] * s.mult.iter.juv    } 
      } 
      
      #  cat(red("wetland is ", wetlands, "\n"))
      # cat(red("popmat.current is ", popmat.current, "\n"))
      #  cat("popmat is ", popmat, "\n")
      
      # instill density dependence for adults  driven by the  number of yr 1s emerging from Berven 2009
      ## This was rejigged to include all adults (noting it will be most heavily influenced by Yr 1s)
      # NOTE K.rel.adult.dem.vec[size] were refined to balance Prmove
      # This is done to offset the very low survival values used
      ## NOTE CHNGED ABOVE FROM >0.5 TO <0.999 CONSIDER THIS IF REJIGGNG ADULTS 
      ##
      totalPop <- sum(n.mat[2:6,i,wetlands])
      K.rel.adult <- (totalPop/K.rel.adult.dem.vec[wetlandMetadata[wetlands,1]]) 
      if (!is.nan(K.rel.adult)) {
        if (K.rel.adult > 0.99)  { 
          K.rel.adult <- 0.99 
        } else if (K.rel.adult <= 0.99 && K.rel.adult > 0.5) {
          for (adgens in 4:6) {
            s.mult.iter.ad <- (as.numeric(coef(fit.expd.adult)[1]) + (K.rel.adult * as.numeric(coef(fit.expd.adult)[2])))^(-1/as.numeric(coef(fit.expd.adult)[3]))
            popmat.current[adgens,(adgens-1)] <- (popmat[adgens,(adgens - 1)] * s.mult.iter.ad)
          } 
        } else if (K.rel.adult <= 0.5)  {
          K.rel.adult <- 0.5 }
      } else {
        stop("The k.rel.adult snippit is misbehaving again")
      }
      #  cat("K.rel.adult is ", K.rel.adult, "\n")
      
      # set popmat.fail for use during dry years
      popmat.fail <- popmat.current
      popmat.fail[1,] <- 0  
  
    
      
      ## Matrix multiplication using the sillHeight for this year to determine Wet or Dry 

           # no ewater so only 
           #fish or no fish
             if (antecedentWet[wetlands] == T && stayWet[1,wetlands] < annualWetness[e,i-1]) {
             popmat.current[1,] <- popmat.current[1,]/600
           }
          #Wet nat
            if (wetlandMetadata[wetlands,2] < sillHeight) {
              antecedentWet[wetlands] <- T
              n.mat[,i+1,wetlands] <- popmat.current %*% n.mat[,i,wetlands] 
           } else {
           #dry nat
             antecedentWet[wetlands] <- F   
             n.mat[,i+1,wetlands] <- popmat.fail %*% n.mat[,i,wetlands]
           }
           
 
      
      
      #  Gen <- Gen + 1
      
      for (ii in 2:6) {       n.mat[ii,i+1,wetlands] <- round(n.mat[ii,i+1,wetlands], digits = 0) }   
      
      #  density feedback for eggs
      #note evidence says there is no density dependence on egg laying in female L. aurea  But there must be some inhibition re maximum volume if nothing else
      # I now use integration. i.e. early in the curve females will lay with 100% success. As it tends towards the pond limit 
      # successive females lay  with diminishing success 
      # above the egg limit for the pond, laying is possible but with a huge inhibition
      K.rel.egg <- (n.mat[ 1,i+1,wetlands]/K.rel.egg.dem.vec[wetlandMetadata[wetlands,1]]) 
      if (!is.nan(K.rel.egg) && (K.rel.egg > 0)) {
        if (K.rel.egg <= 0.99) {  
          intArea <-  integrate(eggfunc,lower=0,upper=K.rel.egg)
          area <- as.numeric(intArea[1])/K.rel.egg
          if (area >= 1) { area <- 1 }
          n.mat[ 1,i+1,wetlands] <- (round(n.mat[1,i+1,wetlands] * area))
        } else if (K.rel.egg > 0.99) {
          # if  > than the limit for the pond then calculate what happens to the first 99% of the pond limit (egg1) then apply the 99th%ile inhibition on the remaining eggs(remain)
          egg1 <- (K.rel.egg.dem.vec[wetlandMetadata[wetlands,1]] * area1)
          remain <- (n.mat[1,i+1,wetlands] - K.rel.egg.dem.vec[wetlandMetadata[wetlands,1]]) 
          n.mat[ 1, i+1, wetlands] <- (round(egg1 + (s.mult.egg.iter * remain)))  
        } else { 
          stop("Crashed at line 1894: the egg conversion value K.rel.egg is misbehaving")
        }  
      }
      
      # theoretically unnecessary but included for stability
      for (clean in 1:6) {
        if (is.nan(n.mat[clean,i+1,wetlands])) {
          n.mat[clean,i+1,wetlands] <- 0 
        }  } 
      
      ## check if I need to increment the patchExt counter
      if (sum(n.mat[,i,wetlands]) > 0 && sum(n.mat[,i+ 1,wetlands]) == 0) {
        patchExt[wetlands] <- patchExt[wetlands] + 1
      }
      
      
      
      
      
      ##  Calculate who moves
      ## remember we are inside the wetland loop. So we will go wetland by wetland and assign to a wetland according to destination
      ## and age
      totalPop <- sum(n.mat[2:6,i + 1,wetlands])
      #  cat("The total population is ", totalPop, "\n")
      emigrants <- 0
      if (totalPop > 0) {
        # we use carrying capacity to determine the probability of moving using Pr(move) = n(t)/K * wetMove 
        # note this is an individual probability of movement, not a proportion of movement
        PrMove <- (totalPop/(K.rel.adult.dem.vec[wetlandMetadata[wetlands,1]])) * wetMove
        ## In order to maintain a PrMove of mean 0.16 +/- 0.14 we constrain the possible values 
        if (PrMove > 0.3) { PrMove <- 0.3 }
        if (PrMove < 0.02) { PrMove <- 0.02 }
        for (potMov in 1:totalPop) { 
          if (runif(1,0,1) < PrMove) { emigrants <- emigrants + 1 }
        }
        #   cat("The Pr.move for this pop is ", PrMove, "\n")
        # cat("The total number of emigrants is  ", emigrants, "\n")         
      }
      
      
      
      if (emigrants >= 1) {
        # Use the determine destinations function to determine destinations        
        destinations <- determine.destination(wetlands,emigrants,adj.wetlandMovement)
        # save this value for error checking
        totalMovers <- length(destinations) - 2 + destinations[length(destinations)] + destinations[length(destinations) - 1]
        
        ## In order to calculate the number of emigrants lost at the upper and lower ends of the reach inside of the determine.destination function
        ## I have hidden them in the returned array. I will now extract them and add them instead to emiLoss NOTE MORTALITY NOT APPLIED
        for (simple in 1:2) {
          emiLoss[1,simple] <- destinations[length(destinations)]
          destinations <- destinations[-(length(destinations))]
        }
        
        ## Do some fancy accounting to assign an age class to each emigrant and remove it from n.mat
        #  cat("a total of ", emigrants, " frogs left from  from wetland ", wetlands, "  NOW we NEED TO ASSIGN AN AGE CLASS \n")
        
        ##calculate the generation division of emigrants and double check that is == emigrants
        ## Can remove this to a function as convenient
        
        ##recreate array 
        emi.gen[] <- 0
        
        #remove from highest to lowest gen, lowest populated gen absorbs the additional emigrants (if required)
        popCheck <- 0
        
        for (sss in 5:1) {
          #  cat("sss goes from 5 to 1. I am at sss = ", sss, "\n")
          
          emi.gen[sss] <- n.mat[sss+1,i+1,wetlands]/totalPop
          emi.gen[sss] <- round(emi.gen[sss] * emigrants, digits = 0)
          if(is.nan(emi.gen[sss])) {
            emi.gen[sss] <- 0
          }
          # doublecheck that the number I am trying to move this generation (emi.gen[sss]) is not larger than the remaining number of emigrants
          if (emi.gen[sss] > (emigrants-popCheck)) {
            emi.gen[sss] <-  (emigrants-popCheck)
          }
          
          popCheck <- popCheck + emi.gen[sss]
          #  cat(blue("popCheck is ", popCheck, "\n"))
          n.mat[sss+1,i+1,wetlands] <- n.mat[sss+1,i+1,wetlands] - emi.gen[sss]
        }
        
        if (popCheck == emigrants) {
          #  cat("the number of emigrants was split appropriately between generations \n")
        } else {
          remainder <- emigrants - popCheck 
          #  cat(red(" the number of emigrants did not split evenly between the generations, I have a remainder of ", (remainder), " attempting resolution \n"))        
          
          for (rem in 2:6) {
            if (n.mat[rem,i+1,wetlands] >= remainder) {
              #   cat(" there are ", n.mat[rem,i+1,wetlands], "frogs in this generation so I will take away ", remainder, "\n")
              n.mat[rem,i+1,wetlands] <-  n.mat[rem,i+1,wetlands] - remainder
              #   cat("which changes it to ", n.mat[rem,i+1,wetlands], " Is this a reduction? I sure hope so!\n")
              emi.gen[rem] <- emi.gen[rem] + remainder
              remainder <- 0
              break
            } 
            if (rem == 6 && remainder > 0) {
              for (remm in 2:6) {
                if (n.mat[remm,i+1,wetlands] > remainder) {   
                  #    cat(red("first run didnt work, I am trialing removing ", remainder, "from n.mat ages 2 through to 6","/n"))
                  n.mat[remm,i+1,wetlands] <-  n.mat[remm,i+1,wetlands] - remainder
                  emi.gen[remm] <- emi.gen[remm] + remainder
                  remainder <- 0
                }
              }
            }
          }
          
          
          # Sorry for this absurd code but hopefully this last loop will catch the errors - when round(emigrants/totalPop) is < 0
          if (sum(emi.gen) != totalMovers) {  
            for (rem in 2:6) {
              if (n.mat[rem,i+1,wetlands] > 0) {
                n.mat[rem,i+1,wetlands] <- n.mat[rem,i+1,wetlands] - 1
                emi.gen[rem-1] <- emi.gen[rem-1] + 1
              }
              if (sum(emi.gen) == totalMovers) { break }
            }
          }
          ## IF THIS IS STILL THROWING ERRORS THEN SHIFT THE ABOVE LOOP TO A FUNCTION AND RECURSE
          
          
          
           # quick error check
          for (rem in 2:6) {
            if (n.mat[rem,i+1,wetlands] < 0 | isFALSE(all.equal(n.mat[rem,i+1,wetlands], as.integer(n.mat[rem,i+1,wetlands])))) {
              stop("Some issues with the removal of emigrants")
            }
          }
        }
        
        ## Current movement survival is a blanket 0.5 Can be changed and moved out of here if necessary
        ## 
        movementSurv <- 0.5
        
        # quick error check 
        #    cat("emi.gen is ", emi.gen, " sum is (", sum(emi.gen), "), and destinations has length ", length(destinations), ". Quick check they are equal . . .\n")
        # cat("totalMovers includes arrivals within the reach and those leaving the reach is NOTE frogs have been removed but not readded", totalMovers, "\n")
        if (sum(emi.gen) != totalMovers) {  stop("length(emi.gen) != length(destinations)")   }
        
        
        ## once emi.gen is complete we use the add.immigrants function to check for survival during movement and then assign to immigrants array
        # we call using immigrants <- add.immigrants(emi.gen,moveSurv,destinations,immigrants) 
        
        ## removed   immigrants <- add.immigrants(emi.gen,moveSurv,destinations,immigrants,successfulMovement,successfulRecolonise,e) 
        
        ## Assign immigrants from emi.gen into immigrants, calculating movement survival on the way past
        # use "immigrants <- add.immigrants(emi.gen,moveSurv,destinations,immigrants) "
        
        
        ## THIS IS NOT THE ACCOUNTING STEP. THIS CHECK SURVIVAL AND ASSIGNS TO THE IMMIGRANTS and successful rescue ARRAYs
        oox <-  length(destinations) 
        #  cat("starting oox loop \n")
        while (oox > 0) {
          #   cat("oox is ", oox, " \n") 
          # check if it survived NOTE ADD STOCHASTICITY TO SURVIVAL FIGURES
          if (runif(1,0,1) <= moveSurv) {
            # cat("survived")
            successfulMovementHold[wetlands,destinations[1]] <- successfulMovementHold[wetlands,destinations[1]] + 1
            ooo <- SelectAge(emi.gen)
            emi.gen[ooo] <- emi.gen[ooo] - 1
            if (sum(n.mat[,i,destinations[oox]]) == 0) {
              successfulRecoloniseHold[wetlands,destinations[oox]] <- successfulRecoloniseHold[wetlands,destinations[oox]] + 1 
            }
            immigrants[destinations[oox],ooo] <- immigrants[destinations[oox],ooo]  + 1
          }
          oox <- oox -1  
          #  cat("finishing oox loop \n")
        }
        ooxrpt <- 0
      } 
      
      
      # call the DetermineAges function to create an array of immigrant ages from outside the reach. Store as immigrantAges. 
      # cat("The sum of emiLoss is ", sum(emiLoss, na.rm = TRUE), "\n")
      if (sum(emiLoss, na.rm = TRUE) > 0) {
        immigrantAges <- DetermineAges(emiLoss)
        #    cat (" Immigrant Ages are ", immigrantAges, "\n")
        if (ncol(immigrantAges) > 0) {
          immigrantDestinations <- immigrantAges
          
          ## assign movement distances
          ## Levy flight use 'dist.rand1 <- sample(dist.vec, 1, replace=T, prob=pr.pred)' 
          for (one in 1:2) {
            for (two in 1:ncol(immigrantDestinations)) {
              if (!is.na(immigrantDestinations[one,two])) {
                dist.rand1 <- sample(dist.vec, 1, replace=T, prob=pr.pred) 
                immigrantDestinations[one,two] <- dist.rand1
              }
            }
          }
          
          # reduce by a distance then assure each one is > 0 then assign a wetland * note mortality during movement has been applied
          for (imi in 1:ncol(immigrantDestinations)) {
            #from DS loop
            stopper <- 2
            dist <- 0
            if (!is.na(immigrantDestinations[1,imi])) {
              dist <- immigrantDestinations[1,imi] - adj.wetlandMovement[1,2]
              if (dist < 0) { dist <- 0 }
              for (immi in 23:2) {
                if (dist < 50) { 
                  immigrantDestinations[1,imi] <- 2 
                  stopper <- 0
                  break
                }
                if (dist < 25) { 
                  immigrantDestinations[1,imi] <- 1 
                  stopper <- 0
                  break
                }
                if (dist > adj.wetlandMovement[1,immi]) {
                  stopper <- immi
                  break
                }
              }
              if (stopper > 0) {
                ifelse(abs(dist - adj.wetlandMovement[1,stopper]) < abs(dist - adj.wetlandMovement[1,stopper+1]), immigrantDestinations[1,imi] <- stopper, immigrantDestinations[1,imi] <- stopper+1)
              }
            }
            
            
            #from US loop
            stopper <- 0
            dist <- 0
            if (!is.na(immigrantDestinations[2,imi])) {
              dist <- abs(immigrantDestinations[2,imi] - adj.wetlandMovement[23,22]) 
              #    cat("Dist1 is ", dist, "\n")
              if (dist < 0) { dist <- 0 }
              #    cat("Dist2 is ", dist, "\n")
              for (immi in 1:22) {
                if (dist <= 150) { 
                  immigrantDestinations[2,imi] <- 23 
                  stopper <- 0
                  break
                }
                if (dist <= 300) {
                  immigrantDestinations[2,imi] <- 22
                  stopper <- 0
                  break
                }
                if (dist > adj.wetlandMovement[23,immi]) {
                  stopper <- immi
                  #    cat("Stopper is ", stopper, "\n")
                  break
                }
              }
              if (stopper > 0) {
                ifelse (abs(dist - adj.wetlandMovement[23,stopper]) < abs(dist - adj.wetlandMovement[23,stopper-1]),immigrantDestinations[2,imi] <- stopper, immigrantDestinations[2,imi] <- stopper-1)
              }
            } 
          }
          ## Store these exogenous immigrants in the exoImmi array for accounting later
          rr <- immigrantAges
          dim(rr) <- c(1,length(rr))
          ss <- immigrantDestinations
          dim(ss) <- c(1,length(ss))
          
          for (steps in 1:length(rr)) {
            if (!is.na(ss[steps]) && !is.na(rr[steps])) {
              ssrr <- c(ss[steps],rr[steps])
              exoImmi <- rbind(exoImmi,ssrr)
              if (max(exoImmi[,1], na.rm = TRUE) > 23) { stop("I have failed to change distance to destination")}
            }
          }
        }  
        # reset the arrays 
        # immigrantAges   
        #  immigrantDestinations
        #  exoImmi
        
        immigrantAges[] <- NA
        immigrantAges <- subset (immigrantAges, select = -c(2:ncol(immigrantAges)))    
        immigrantDestinations[] <- NA
        #  immigrantDestinations <- subset (immigrantDestinations, select = -c(2:ncol(immigrantDestinations)))    
        emiLoss[] <- 0 
        
      }
      
      totalPop <- sum(n.mat[2:6,i + 1,wetlands])
      #   cat("The total population is ", totalPop, "\n")
      
      
      ##  Last line of the 'wetlands' which cycles through each of our wetlands in turn
    }
    
    
    ##  not sure if this will work here or if it needs to move up a loop
    ##  add all the immigrants to the respective n.mats and then reset the immigrants array
    if (sum(immigrants > 0)) {  
      #  cat(blue("I am adding the new immigrants now \n")) 
      for (kk in 1:wetlandNum) {
        for (kkk in 1:5) {
          if (immigrants[kk,kkk] > 0) {
            #  cat(immigrants[kk,kkk], " frogs which were age class ", kkk+1, " successfully arrived at wetland ", kk, "changing the population from ",  n.mat[kkk+1,i+1,kk], " to ", (n.mat[kkk+1,i+1,kk] + immigrants[kk,kkk]), "\n")  
            n.mat[kkk+1,i+1,kk] <- (n.mat[kkk+1,i+1,kk] + immigrants[kk,kkk])
          }
        }  
      }
    } else {
      #     cat(blue (" no immigrants this year \n"))
    }
    #   print("Step3 Complete \n")
    ## RESET the temporary arrays 
    
    immigrants[] <- 0 
    #   adj.wetlandMovement[] <- NA
    
    ## Now add in the exogenous immigrants but don't track recolonisations     
    if (sum(exoImmi, na.rm = TRUE) > 0) {
      for (trax in 2:nrow(exoImmi)) {
        n.mat[exoImmi[trax,2]+1,i+1,exoImmi[trax,1]] <- n.mat[exoImmi[trax,2]+1,i+1,exoImmi[trax,1]] + 1
      }
      exoImmi <- exoImmi[1,]  
    }
    
    #  cat(" i is ", i, "\n")
    ## On the last generation tally up all the population totals and adult totals and ewaterevents
  #  if (i == (generations)) {
   #   for (cycle in 1:wetlandNum) {
    #    finalControl[cycle,e] <- sum(n.mat[,i,cycle])
     #   finalControl[cycle+23,e] <- sum(n.mat[2:6,i,cycle])
        # THE EWATER NUMBERS GO HERE
    #  }
  #  }
    
    ## and a quick error check for fun
    if (i > (generations)) {
      stop("failed to break at the 23rd year")
    }
 
    adultHolder <- 0
    ## THE ACCOUNTING LOOP TO STORE THE VALUES
    if (i %% 10 ==0) {
      colmn <- (i/10) 
    for (l in 1:23) {
      #store allPop wetland
      allPop[l,colmn + 1] <- allPop[l,colmn + 1] + sum(n.mat[1:6,i,l])
      
    #add value to row 24 (total pops alive)  
    if (sum(n.mat[1:6,i,l] > 0)) { 
      allPop[24,colmn + 1]  <- allPop[24,colmn + 1] + 1  
      occupancyHeatC[l,colmn + 1,SDTracker] <- occupancyHeatC[l,colmn + 1,SDTracker] + 1 
      }
      
   #store adultPop by wetland
   adultPop[l,colmn + 1] <- adultPop[l,colmn + 1] + sum(n.mat[2:6,i,l])
   adultHolder <- adultHolder + sum(n.mat[2:6,i,l])                                                           
      #add value to row 24 (adult pops alive)
      if (sum(n.mat[2:6,i,l] > 0)) { adultPop[24,colmn + 1]  <- adultPop[24,colmn + 1] + 1  }
    }  
      
      SDArray[1,colmn,SDTracker] <- adultHolder
      adultHolder <- 0
    
       # store eWater event by wetland
      eWater[1:23,colmn] <- eWater[1:23,colmn] + eWaterHold
      SDArray[2,colmn,SDTracker] <- sum(eWaterHold)
      
      # store successful movements 23 * 23 (from (row) to (column)
      successfulMovement[,,colmn] <- successfulMovement[,,colmn] + successfulMovementHold[,]
      SDArray[4,colmn,SDTracker] <- sum(successfulMovementHold[,])
      
      # store successful recolonisations 23 * 23 (from (row) to (column)
      successfulRecolonise[,,colmn] <- successfulRecolonise[,,colmn] + successfulRecoloniseHold[,]
      SDArray[5,colmn,SDTracker] <- sum(successfulRecoloniseHold[,])
      
      patchExtAll[1:23,colmn + 1] <- patchExtAll[1:23,colmn + 1] + patchExt[]
      SDArray[3,colmn,SDTracker] <- sum(patchExt[])
      
      for (l in 1:23) {
        SDWetlandMove[l,colmn,1,SDTracker] <- sum(successfulMovementHold[l,])
        SDWetlandMove[l,colmn,2,SDTracker] <- sum(successfulRecoloniseHold[l,])
      }
      
     
    }  
   
  
    
    
       
    
    ## Last line of the generation loop which controls how many years we run for t is (generations + 10)
  }
#  for (rroww in 1:23) {
#    successfulRecolonise[rroww,24,e] <- sum(successfulRecolonise[rroww,1:23,e])
#  }
  ## Last line of Outermost Loop which is the Iteration Loop (currently set to one)
  
  
}


write.csv(occupancyHeatC,"C:/Workspace/metaOutput/ContOccupancyHeat.csv",row.names = FALSE)










library(ggplot2)
library(plotly) 

## REFORMAT THE DATA FOR A RASTER MAP
aaa <- array(data = 0, dim = (23*7))
dim(aaa) <- c(23,7)
for (I in 1:7) {
  for (II in 1:23) {
    aaa[II,I] <- sum(occupancyHeatSSTH[II,I,1:iter])/iter
  }
}
AAA <- as.data.frame(aaa)

colnames(AAA) <- c(0,10,20,30,40,50,60)
WL <- rep(c(1:23),times = 7)
# Transform data to long form
aaa.melt <- melt(AAA)
BBB <- as.data.frame(aaa.melt,WL)
BBB$wetlands=WL
colnames(BBB) <- c("years", "PrSv","wetlands")
#BBB <- as.numeric(BBB)
BBB$years <- as.numeric(BBB$years)
BBB$years <- BBB$years*10
BBB$years <- BBB$years-10
BBB$PrSv[1:23] <- startPrSv

## CALCULATE THE MN & SD FOR A COMPARATIVE LINE GRAPH
VVV <- c(1:iter)
VVV[] <- 0

# treatment 0 - control
CPopSvMn <- c(0,0,0,0,0,0,0)
CPopSvUp <- c(0,0,0,0,0,0,0)
CPopSvDwn <- c(0,0,0,0,0,0,0)
for (j in 1:6) {
  for (i in 1:iter) {
    VVV[i] <- sum(occupancyHeatC[,j+1,i])
  }
  CPopSvMn[j + 1] <- mean(VVV)
  CPopSvUp[j + 1] <- CPopSvMn[j + 1] + SD(VVV)
  CPopSvDwn[j + 1] <- CPopSvMn[j + 1] - SD(VVV)
}

# treatment 1 SSTH
SSTHPopSvMn <- c(0,0,0,0,0,0,0)
SSTHPopSvUp <- c(0,0,0,0,0,0,0)
SSTHPopSvDwn <- c(0,0,0,0,0,0,0)
for (j in 1:6) {
for (i in 1:iter) {
  VVV[i] <- sum(occupancyHeatSSTH[,j+1,i])
}
  SSTHPopSvMn[j + 1] <- mean(VVV)
  SSTHPopSvUp[j + 1] <- SSTHPopSvMn[j + 1] + SD(VVV)
  SSTHPopSvDwn[j + 1] <- SSTHPopSvMn[j + 1] - SD(VVV)
  }

axs <- c(0,0,0,0,0,0,0)
years <- c(0,10,20,30,40,50,60)
SSTHPopSvdf <- data.frame(years,axs,SSTHPopSvMn,SSTHPopSvUp,SSTHPopSvDwn,CPopSvMn,CPopSvUp,CPopSvDwn) 


SSTHPopSvTEST <- ggplot(data=SSTHPopSvdf, aes(x=years,y=), colour=Flow) + 
  # ggplot(data=NatExt) +
  theme(axis.line=element_line(color="black",size=1)) +
  theme_classic(base_size = 22) +
  # scale_x_discrete(breaks = 1:9,labels = Heights) +
  scale_x_continuous(limits = c(0, 60),breaks=seq(0,60,10),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 23),breaks=seq(0,23,5), expand = c(0, 0)) +
  # draw_image(Piccy, x = 8, y = -19.23, width = 40, height = 40, scale = 1,clip = "inherit", interpolate = TRUE, hjust = 0, vjust = 0) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank()) + 
  ggtitle("SvPops +/-SD)") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  #labs(y= "adult Pop Increase", x="years") + 
  theme(aspect.ratio=.7) +
  theme(axis.text.y = element_text(size=12, hjust= -0.5)) +
  
  #  geom_ribbon(aes(ymin = Dwn2, ymax = Up2), fill = "grey90") +
  #geom_ribbon(aes(ymin = PercAdDwn1, ymax = PercAdUp1), fill = "grey90") +
  geom_ribbon(aes(ymin = SSTHPopSvDwn, ymax = SSTHPopSvUp), fill = "grey90") +
    geom_ribbon(aes(ymin = CPopSvDwn, ymax = CPopSvUp), fill = "grey90") +
  #  geom_ribbon(aes(ymin = Dwn4, ymax = Up4), fill = "grey90") +
 # geom_line(linetype="solid",color="black",size=1,aes(x=years,y=axs), show.legend = TRUE) +
  geom_line(linetype="solid",color="blue",size=1,aes(x=years,y=CPopSvMn), show.legend = TRUE) +
  geom_line(linetype="solid",color="black",size=1,aes(x=years,y=SSTHPopSvMn), show.legend = TRUE) 
#  geom_line(linetype="longdash",color="black",size=1,aes(x=years,y=PercAdMn2), show.legend = TRUE) +
 # geom_line(linetype="solid",color="red",size=1,aes(x=years,y=PercAdMn4), show.legend = TRUE) +
  # geom_line(linetype="dotted",color="black",size=1.3,aes(x=years,y=PercAdMn3), show.legend = TRUE) 


SSTHPopSvTEST


#  stat_contour()

(ggplot(data=BBB, aes(x=years, y=wetlands, z=PrSv))
 # + geom_point(aes(colour=PrSv))
#  + stat_density2d()
#+  geom_density2d_filled(aes(colour=PrSv))
#  + geom_bin2d(binwidth = c(0.1,1))
  + geom_raster(aes(fill=PrSv))
 # + geom_tile(aes(fill=PrSv))
 # +  geom_hex(aes(colour=PrSv))
)

## calculate PrSv at year0

holder <- c(1:10000)
holder[] <- 0
startPrSv <- 1:23
startPrSv[] <- 0
for (j in 1:23) {
for (i in 1:10000) {
    r <- sum(startadPopsAll[(j*6):(j*6)-5,i])
   if (r > 0) { holder[i] <- holder[i] + 1 } 
  }
  startPrSv[j] <- sum(holder)/10000 
  holder[] <- 0
}




beep(3)
## finalControl is the output.csv for this run
## 23*3 rows  23*1 is totalPop, 23*2 is adultPop, 23*3 is eWaterInterventions
write.csv(successfulRecolonise,"C:/Workspace/metaOutput/ContRecol.csv",row.names = FALSE)
write.csv(successfulMovement,"C:/Workspace/metaOutput/ContMovement.csv",row.names = FALSE)
write.csv(eWater,"C:/Workspace/metaOutput/ConteWater.csv",row.names = FALSE)
write.csv(adultPop,"C:/Workspace/metaOutput/ContadultPop.csv",row.names = FALSE)
write.csv(allPop,"C:/Workspace/metaOutput/ContallPop.csv",row.names = FALSE)
write.csv(patchExtAll,"C:/Workspace/metaOutput/ContpatchExt.csv",row.names = FALSE)
write.csv(SDArray,"C:/Workspace/metaOutput/ContSDArray.csv",row.names = FALSE)
write.csv(SDWetlandMove,"C:/Workspace/metaOutput/ContSDmovement.csv",row.names = FALSE)

write.csv(occupancyHeatC,"C:/Workspace/metaOutput/ContOccupancyHeat.csv",row.names = FALSE)




