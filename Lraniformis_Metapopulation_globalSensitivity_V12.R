# The first line ####
##############################################################################################################
##############################################################################################################
## Short run attenuated code designed for sensitivity analysis using a Latin Hypercube
## Uses Large wetlands under regulated flow conditions at a sill height of 7.5 m
##############################################################################################################
# Nesting is 1000 runs then wet between the dry then duration of the dry then population size then T years  

## Variables for sensitivity
# clutch.size.lr <- 1000 to 6000
# hatch.dur <- 0.5 to 14
# hatch.pr <- 0.5 - .98
# tadpole.dur <- 50 to 90
# tadpole.mn <- 0.03 to 0.6
# ad.s.yr.mn <- 0.03 <- 0.35
# maxAge <- 5 to 10
# K.rel.tad.dem.vec[3] <- 500 to 5000
# K.rel.juv.dem.vec[3] <- 200 to 5000

## CJA Bradshaw & Rupert Mathwin
## June 2021

## Remove everything
rm(list = ls())

PrWL0 <- c(1:23)

##Changed source file location below
source("C:/workspace/math0286/R/win-library/3.6/matrixOperators.r")

## libraries
library(adehabitatLT)
library(beepr)
library(crayon)
library(data.table)
library(DataCombine)
library(DescTools)
library(dismo)
library(doSNOW)
library(foreach)
library(gbm)
library(grid) #, lib.loc = "C:/Program Files/R/R-4.0.2/library")
library(iterators)
library(lhs)
library(parallel)
library(plyr)
library(pracma)
library(raster)
library(readr)
library(reshape2)
library(snow)
library(VGAM)

# setwd("~/Documents/Papers/Amphibians/bell frog pop model/analysis/gsa")

## source (update when appropriate)
# source("matrixOperators.r") 


## Set up parallel processing (nproc is the number of processing cores to use)
cores <- detectCores()
nproc <- (cores - 2)
cl.tmp = makeCluster(rep('localhost', nproc), type='SOCK')
registerDoSNOW(cl.tmp)
getDoParWorkers()

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


## encapsulate the core metapopulation model as a function for parallelisation
sbf_sim <- function(input, dir.nm, rowNum) {
  
  ## assign all parameter values
  for (d in 1:ncol(input)) {assign(names(input)[d], input[,d])}
   # complementary log-log
  cloglog <- function(x) log(-log(1-x))
  
    ####################################################
    ## necessary input calculations for stochastic model
    ####################################################
   
  antecedentWet <- logical(length=23)
  iter <- 10
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
  
  #############################################################################################################
  ## Custom Functions Live Below
  #############################################################################################################
  # beta distribution shape parameter estimator function
  ##Generates an alpha and a beta value to inform the beta distribution 
  estBetaParams <- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta)) }
  
  ## Levy flight Setup code
  ## use 'dist.rand1 <- sample(dist.vec, 100, replace=T, prob=pr.pred)' to produce 100 values
  ## to create dist.rand1
  # also the plot window has to be a certain size so if it crashes expand the Rstudio plot window
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
     # quick reminder c is US loss, d is DS loss
    moveDist <- as.numeric(moveDist[!is.na(moveDist)])
    moveDist <- append(moveDist,c)
    moveDist <- append(moveDist,d)
    return(moveDist)
  }
  
  
  ## Assign immigrants from emi.gen into immigrants, calculating movement survival on the way past
  # call using "immigrants <- add.immigrants(emi.gen,moveSurv,destinations,immigrants) "
  add.immigrants <- function(emigrntGen,survival,desti,immiArray,storeMove,storeRecol,it) { 
    for (oo in 1:length(emigrntGen)) {
      while (emigrntGen[oo] > 0) {
        emigrntGen[oo] <- emigrntGen[oo] - 1
        # check if it survived
        if (runif(1,0,1) <= survival) {
          immiArray[desti[1],oo] <- immiArray[desti[1],oo] + 1
          desti <- desti[-1]
        }
      }
    }
    return(immiArray)
  }
  
  
  ## determine age of immigrants DetermineAges(emiLoss)
  #code snippet to assign ageClass according to expected age structure (noting range 1 is 1-2yrs, 2 is 2-3 yrs)
  DetermineAges <- function(immiNums) {
    fromUS <- NA
    fromDS <- NA
    counter <- 1 
    for (sec in 1:2) {
      while (counter <= immiNums[1,sec]) {
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
    return(immiAges)
  }
  
  # this is constrained later but retained here
  ## set time limit for projection in 1-yr increments
  yr.now <- 2020
  #************************
  #yr.end <- 2020 + round((40*gen.l), 0) # set projection end date
  yr.end <- 2020 + 86
  #************************
  ##Note I have constrained t to 60 years for this model
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
 tomet.dur.mn <- round(sum(c(mean(hatch.dur), mean(tadpole.dur))), 0)
  toad.dur.mn <- 365 - tomet.dur.mn
   # annual adult survival from Pickett et al (L aurea)
  ad.s.yr.mn <- 0.2172 # Litoria aurea (Pickett et al.)
  ad.s.yr.sd <- 0.087 # Litoria aurea (Pickett et al.)
  ad.s.yr.alpha <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$alpha
  ad.s.yr.beta <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$beta
  
  ##Calculate our survivals
  #likelyhood of egg hatching and tadpole surviving to metamorphosis
  tomet.s.iter <- (rbeta(1, tp.s.alpha, tp.s.beta)) * (runif(1, min=0.933, max=1))
  #sample a daily probability of survival
  toad.daily.s.iter <- nthroot(rbeta(1,ad.s.yr.alpha, ad.s.yr.beta) , 365)
  toad.s.season.iter <- toad.daily.s.iter ^ toad.dur.mn
  toad.s.iter <- tomet.s.iter * toad.s.season.iter                     
  #Create the temp survival vector, to adult survival then 4 adult survivals
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
  DD.dat <- data.frame(K.egg.vec, surv.mult.egg.vec)
  param.init <- c(1.01, 9, 2)
  fit.expd.egg <- nls(surv.mult.egg.vec ~ (a - ((K.egg.vec ^ (b/K.egg.vec)) * (K.egg.vec ^ d))),
                      algorithm = "port",
                      start = c(a = param.init[1], b = param.init[2], d = param.init[3]),
                      trace = TRUE,      
                      nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
  K.pred.egg.vec <- seq(0.01, 1, 0.01)
  pred.surv.egg.mult <- (as.numeric(coef(fit.expd.egg)[1]) - (K.pred.egg.vec ^ (as.numeric(coef(fit.expd.egg)[2])/K.pred.egg.vec) * (K.pred.egg.vec ^ (coef(fit.expd.egg)[3]))))
  eggfunc <- function(x) (1.01 - ((x ^ (9/x)) * (x ^ 2)))
  intArea <- integrate(eggfunc,lower=0,upper=0.99)
  area1 <- (as.numeric(intArea[1])/0.99)
  s.mult.egg.iter <- (as.numeric(coef(fit.expd.egg)[1]) - (0.999 ^ (as.numeric(coef(fit.expd.egg)[2])/0.999) * (0.999 ^ (coef(fit.expd.egg)[3]))))
  
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
  K.pred.tad.vec <- seq(K.lo,1,0.01)
  pred.surv.tad.mult <- (as.numeric(coef(fit.expd.tad)[1]) + (K.pred.tad.vec * as.numeric(coef(fit.expd.tad)[2])))^(-1/as.numeric(coef(fit.expd.tad)[3]))
  
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
  K.pred.juv.vec <- seq(K.lo,1,0.01)
  pred.surv.juv.mult <- (as.numeric(coef(fit.expd.juv)[1]) + (K.pred.juv.vec * as.numeric(coef(fit.expd.juv)[2])))^(-1/as.numeric(coef(fit.expd.juv)[3]))
  
  fit.expd.adult <- fit.expd.tad
  
  ### ###
  ## Populate underlying movement data (site-to-site distance which incorporates site-to-site movement type)
  ### ###
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
  
  ### Create an array of distance between wetlands where the distance is adjusted to incorporate landscape resistance downstream easier/shorter, upstream harder/further.
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
  
  ## Create a holder for emigration from the reach (top and bottom) during dispersal
  emiLoss <- array(data = 0, dim = 2)
  dim(emiLoss) <- c(1,2)
  colnames(emiLoss) <- c("DS Loss","US Loss")
  
  exoImmi <- array(data = NA, dim = 2)
  dim(exoImmi) <- c(1,2)
  colnames(exoImmi) <- c("Wetland","AgeClass")
  
  ## AdultPop and AllPop will have an additional row 
  ## Row24 is the cumulative total of how many pops were alive at each 10yr step
  ## Create an array to hold all of the final populations at end of each iter/10
  adultPop <- array(data = 0, dim = 24 * 7)
  dim(adultPop) <- c(24,7)
  
  allPop <- array(data = 0, dim = 24 * 7)
  dim(allPop) <- c(24,7)
  
   # SDWetlandMove will track the number of movements and recolonisations from each wetland 
  # dims are 23wetlands,7 time steps,2 (1move,2recol),10000 iters
  SDWetlandMove <- array(data = 0, dim = 23*7*2 * 10000)
  dim(SDWetlandMove) <- c(23,7,2,10000)
  
  wetCycle <- array(data=0,dim=23*5)
  dim(wetCycle) <- c(5,23)
  
  # Control
  eWaterSites <- c(0)
  
  banrockFish <- c(F,T,T)
  
  ## KEEP THIS   Array to hold for probability of occupancy ?heatmap
  occupancyHeatC <- array(data=0,dim=23*7*10000)
  dim(occupancyHeatC) <- c(23,7,10000)
  
    ###############################################################################################################
    ## And NOW for the actual model!
    ###############################################################################################################
  
  print("two --")

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
  
   SDTracker <- 0
  immigrantDestinations <- 0
  
  ## The Outermost Loop: iterate the process iter times (spoiler alert it's one time)
  ## note I currently have 150 unique centuries if it exceeds 150 we will have to start again
  for (e in 1:iter) {
    cat("starting iteration ", e, "\n")  
    SDTracker <- SDTracker + 1
    sillCounter <- 0 
    antecedentWet[] <- F
    wetCycle[] <- 0
    
    ## The Second Loop: run the current projection set up for the number of generations (years actually) + 10 for burn-in 
    ## NOTE I have generated centuries so this value can go up to 90 yrs as required, currently at 50
    for (i in 1:(generations)) {
      ## work out the wetness figures for this year the sillHeight and wetMod for this year
      # note I have applied the reduction of 10cm to allow for the threshold in this step
      sillCounter <- sillCounter + 1
      sillHeight <- sillForecast[e,sillCounter]
    
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
        }
      }
     
      #iterate the runCounter used for errorchecking
      if (i==1) { runCounter[wetlandNum] <- (runCounter[wetlandNum] + 1) 
      }
      
      # if there are no frogs alive break from this loop
      if (sum(n.mat[,i,]) == 0) {
        break  
      }
    
        ## Cycle through each of the wetlands, The Wetlands Loop
      for (wetlands in 1:wetlandNum) {
       
        # resample the duration of egg to tadpole
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
        
        # instil density dependence for tadpoles growing into year 1 adults
        if (wetlandMetadata[wetlands,1] == 2) {
         K.rel.tad.use <-  K.rel.tad.Hold2 
        } else if (wetlandMetadata[wetlands,1] == 3) {
          K.rel.tad.use <-  K.rel.tad.Hold3 
        } else if (wetlandMetadata[wetlands,1] == 4) {
          K.rel.tad.use <-  K.rel.tad.Hold4 
        } else {
          stop("line 965 K.rel.tad.Hold")
        }
        
        K.rel.tad <- (n.mat[1,i,wetlands]/K.rel.tad.use)
        if (!is.nan(K.rel.tad)) {
          if (K.rel.tad > 2.1)  { K.rel.tad <- 2.1 }
          if (K.rel.tad <= 2.1) {
            s.mult.iter.tad <- (as.numeric(coef(fit.expd.tad)[1]) + (K.rel.tad * as.numeric(coef(fit.expd.tad)[2])))^(-1/as.numeric(coef(fit.expd.tad)[3]))
            popmat.current[2,1] <- popmat[2,1] * s.mult.iter.tad   } 
        } 
        
        
        if (wetlandMetadata[wetlands,1] == 2) {
          K.rel.juv.use <-  K.rel.juv.Hold2 
        } else if (wetlandMetadata[wetlands,1] == 3) {
          K.rel.juv.use <-  K.rel.juv.Hold3 
        } else if (wetlandMetadata[wetlands,1] == 4) {
          K.rel.juv.use <-  K.rel.juv.Hold4 
        } else {
          stop("line 965 K.rel.tad.Hold")
        } 
        
        # instil density dependence for juveniles (1 - 2 years) is driven by the  number of yr 1 present/competing per 
        K.rel.juv <- (n.mat[2,i,wetlands]/K.rel.juv.use) 
        if (!is.nan(K.rel.juv)) {
          if (K.rel.juv > 2.1)  { K.rel.juv <- 2.1 }
          if (K.rel.juv <= 2.1) {
            s.mult.iter.juv <- (as.numeric(coef(fit.expd.juv)[1]) + (K.rel.juv * as.numeric(coef(fit.expd.juv)[2])))^(-1/as.numeric(coef(fit.expd.juv)[3]))
            popmat.current[3,2] <- popmat[3,2] * s.mult.iter.juv    } 
        } 
        
        # instill density dependence for adults  driven by the  number of yr 1s emerging from Berven 2009
        
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
        
        ##  Calculate who disperses
        ## remember we are inside the wetland loop. So we will go wetland by wetland and assign to a wetland according to destination (retains age)
        totalPop <- sum(n.mat[2:6,i + 1,wetlands])
        emigrants <- 0
        if (totalPop > 0) {
          # This version we assign PrMove as the Latin Hypercube and not using carrying capacity
          # note this is an individual probability of movement, not a proportion of movement
          for (potMov in 1:totalPop) { 
            if (runif(1,0,1) < PrMove) { emigrants <- emigrants + 1 }
          }
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
          
          ##recreate array 
          emi.gen[] <- 0
          
          #remove from highest to lowest gen, lowest populated gen absorbs the additional emigrants (if required)
          popCheck <- 0
          
          for (sss in 5:1) {
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
            n.mat[sss+1,i+1,wetlands] <- n.mat[sss+1,i+1,wetlands] - emi.gen[sss]
          }
          
          if (popCheck == emigrants) {
          } else {
            remainder <- emigrants - popCheck 
            for (rem in 2:6) {
              if (n.mat[rem,i+1,wetlands] >= remainder) {
                n.mat[rem,i+1,wetlands] <-  n.mat[rem,i+1,wetlands] - remainder
                emi.gen[rem] <- emi.gen[rem] + remainder
                remainder <- 0
                break
              } 
              if (rem == 6 && remainder > 0) {
                for (remm in 2:6) {
                  if (n.mat[remm,i+1,wetlands] > remainder) {   
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
            
            # quick error check
            for (rem in 2:6) {
              if (n.mat[rem,i+1,wetlands] < 0 | isFALSE(all.equal(n.mat[rem,i+1,wetlands], as.integer(n.mat[rem,i+1,wetlands])))) {
                stop("Some issues with the removal of emigrants")
              }
            }
          }
          
          movementSurv <- 0.5
          
          # quick error check 
          if (sum(emi.gen) != totalMovers) {  stop("length(emi.gen) != length(destinations)")   }
          
          ## THIS IS NOT THE ACCOUNTING STEP. THIS CHECK SURVIVAL AND ASSIGNS TO THE IMMIGRANTS and successful rescue ARRAYs
          oox <-  length(destinations) 
          while (oox > 0) {
            # check if it survived NOTE ADD STOCHASTICITY TO SURVIVAL FIGURES
            if (runif(1,0,1) <= moveSurv) {
              ooo <- SelectAge(emi.gen)
              emi.gen[ooo] <- emi.gen[ooo] - 1
              
              immigrants[destinations[oox],ooo] <- immigrants[destinations[oox],ooo]  + 1
            }
            oox <- oox -1  
          }
          ooxrpt <- 0
        } 
        
        
        # call the DetermineAges function to create an array of immigrant ages from outside the reach. Store as immigrantAges. 
        if (sum(emiLoss, na.rm = TRUE) > 0) {
          immigrantAges <- DetermineAges(emiLoss)
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
                if (dist < 0) { dist <- 0 }
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
          emiLoss[] <- 0 
          
        }
        totalPop <- sum(n.mat[2:6,i + 1,wetlands])
        
        ##  Last line of the 'wetlands' which cycles through each of our wetlands in turn
      }
      
      ##  add all the immigrants to the respective n.mats and then reset the immigrants array
      if (sum(immigrants > 0)) {  
        for (kk in 1:wetlandNum) {
          for (kkk in 1:5) {
            if (immigrants[kk,kkk] > 0) {
              n.mat[kkk+1,i+1,kk] <- (n.mat[kkk+1,i+1,kk] + immigrants[kk,kkk])
            }
          }  
        }
      } else {
      }
      
      ## RESET the temporary arrays 
      immigrants[] <- 0 
      
      ## Now add in the exogenous immigrants but don't track recolonisations     
      if (sum(exoImmi, na.rm = TRUE) > 0) {
        for (trax in 2:nrow(exoImmi)) {
          n.mat[exoImmi[trax,2]+1,i+1,exoImmi[trax,1]] <- n.mat[exoImmi[trax,2]+1,i+1,exoImmi[trax,1]] + 1
        }
        exoImmi <- exoImmi[1,]  
      }
      
      if (i == generations) {
        for (j in 1:23) { 
          #Store this for later
          if (n.mat[1:6,i,j] > 0) {
            occupancyHeatC[j,7,1] <- occupancyHeatC[j,7,1] + 1
      }
        }
      }
       
       ## and a quick error check for fun
      if (i > (generations)) {
        stop("failed to break at the 'generations' year")
      }
           
           ##  Last line of the generation loop (86 years)
    }
    
          ##  Last line of the 'iter' Loop  
  }
  
  PrSvReach <- sum(occupancyHeatC[1:23,7,1])/iter/23
  
  # save
  input$PrExt <- PrSvReach
  save.nm <- paste0('res',sprintf("%09.0f", rowNum))
  assign(save.nm, input)
  save(list=save.nm,file=paste(dir.nm,save.nm,sep='/'))
  
    print("*******************")
    print(d) 
    print("*******************")
} # end d lh loop

## parameter ranges
ranges <- list()

ranges$clutch.size.lr <- c(1000,6000)
ranges$hatch.dur <- c(0.5, 5)
ranges$hatch.pr <- c(0.5, 0.98)
ranges$tadpole.dur <- c(50, 90)
ranges$tadpole.mn <- c(0.03, 0.6)
ranges$ad.s.yr.mn <- c(0.03, 0.35)
ranges$K.rel.tad.Hold2 <- c(150, 1000)
ranges$K.rel.tad.Hold3 <- c(400, 2000)
ranges$K.rel.tad.Hold4 <- c(500, 3000)
ranges$K.rel.juv.Hold2 <- c(100, 700)
ranges$K.rel.juv.Hold3 <- c(300, 2000)
ranges$K.rel.juv.Hold4 <- c(500, 3000)
ranges$moveSurv <- c(0.1, 0.7)
ranges$PrMove <- c(0.05, 0.4)
ranges$max.dist <- c(5000,50000)

## create hypercube+
nSamples <- 10000
lh <- data.frame(randomLHS(n=nSamples, k=length(ranges)))
names(lh) <- names(ranges)

## convert parameters to required scale
for (j in 1:ncol(lh)) {
  par <- names(lh)[j]
  lh[,par] <- qunif(lh[,j], min=ranges[[par]][1], max=ranges[[par]][2]) ## continuous
}

## number of iterations for each parameter set
lh$iter <- 10

## folder for saving the results of each row
## we could just store in memory, but then if something breaks we will lose the lot
dir.nm <- 'MetaPopGSAsbfTestFolder'
dir.create(dir.nm)

## run in parallel
res <- foreach(rowNum=1:nrow(lh),.verbose=T) %do% {sbf_sim(input=lh[rowNum,],dir.nm=dir.nm,rowNum=rowNum)}

## retrieve results
res.nms <- list.files(dir.nm)
res.list <- lapply(res.nms, function(x) {load(paste(dir.nm,x,sep='/')) ; print(x) ; return(eval(as.name(x)))})
#res.list1 <- lapply(res.nms1, function(x) {load(paste('GSA50kSparatest5Dmax50pc1',x,sep='/')) ; print(x) ; return(eval(as.name(x)))})
dat <- data.frame(rbindlist(res.list))
#dat <- rbind(dat1,dat2,dat3,dat4)
head(dat)
dim(dat)[1]
tail(dat)
sum(is.na(dat$PrExt))



#########
## BRT ####
#########
dat.nona <- data.frame(na.omit(dat[!is.infinite(rowSums(dat)),]))
dat.nona <- dat.nona[,-9]
dat.nona$cllPrExt <- cloglog(dat.nona$PrExt)
dim(dat.nona)[1]
head(dat.nona)




brt.fit <- gbm.step(dat.nona, gbm.x = attr(dat.nona, "names")[1:14], gbm.y = attr(dat.nona, "names")[16], family="gaussian", max.trees=100000, 
                    tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.75, tree.complexity = 2)

summary(brt.fit)
dim(dat.nona)[1]
D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)

CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
CV.cor
CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
CV.cor.se
print(c(CV.cor, CV.cor.se))

eq.sp.points <- 100
RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=14)
## output average predictions
for (p in 1:8) {
  RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
  RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
}
RESP.val.dat <- as.data.frame(RESP.val)
colnames(RESP.val.dat) <- brt.fit$var.names
RESP.pred.dat <- as.data.frame(RESP.pred)
colnames(RESP.pred.dat) <- brt.fit$var.names
RESP.val.dat
RESP.pred.dat

## plot highest-influence variables only
plot(RESP.val.dat[,6], RESP.pred.dat[,6], type="l", xlab="mean adult annual survival", ylab="Pr(ext)")
plot(RESP.val.dat[,5], RESP.pred.dat[,5], type="l", xlab="tadpole survival to metamorphosis", ylab="Pr(ext)")

setwd("C:/workspace/math0286/R/")
write.table(RESP.val.dat,file="BRT.val.MetaGSAsbf.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.pred.dat,file="BRT.pred.MetaGSAsbf.csv",sep=",", row.names = T, col.names = T)

save.image(paste("MetaGSAsbf",".RData",sep=""))
