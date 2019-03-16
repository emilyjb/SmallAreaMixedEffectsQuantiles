rm(list=ls(all=TRUE))

### Set your working directory to the folder that contains the files posted to Github. For example:
setwd("C:/Users/Emily/Documents/GitHub/SmallAreaMixedEffectsQuantiles")


###  Load libraries
library(nlme)
library("lqmm")
library("sae")
library("survey")
library("quantreg")
library("sn")
library("eva")
 
### Source the following R files
source("storoutputsimpleEMsim.R") 		# Store output
source("genpopmixedllcomb.R")			# Generate data
source("NEBfuns.R")				# For EBP predictors
source("LIGPDfuns.R")				# For LIGPD predictors
source("Bootfuns.R")				# For Bootstrap MSE computation

### Population and sample size configuration:
N <- 20000
nis <- rep(c(rep(5,5),rep(10,5), rep(20,5)),times = 4)
Nis <- round(N*nis/sum(nis),0)
CNis <- cumsum(Nis)
D <- length(nis)
n <- sum(nis)
N <- sum(Nis)

### Bootstrap sample size (denoted by T in the manuscript)
B <- 100

# Sampling fractions
pis  <- nis/Nis
pissamp <- rep(nis/Nis, times = nis)
Nispop <- rep(Nis, Nis)
areafac.pop <- rep(1:D, Nis)
Nissamp <- rep(Nis, nis)

GN <- model.matrix(lm(Nispop~as.factor(areafac.pop)-1))

### Set parameters:

beta0 <- -1.5; beta1 <- 0.5
sig2le <- 1 ; sig2lu <- 0.5
mulx <- 0
sig2lx <- 1

####  Set distribution options:
edistfun <- SNgenfun ##( Other option used in manuscript is "chigenfunH" used for chi-sqaure)
e.parms <- -5 # For chi-square used in manuscript, change e-distribution parameters to (2, 0.1) 
bdistfun <- laplacegenfun ##(Options in manuscript are "laplacegenfun" for "Laplace" and "normalbgenfun" for "Normal")
b.dist <- "Laplace" ##(Change to "Normal" for normal bi.)
doBoot <- TRUE  ##(Change to "FALSE" to skip running bootstrap.)

time.start.all <- Sys.time()

options("contrasts" = c("contr.sum", "contr.poly"))  ### Use sum to zero constraints for initial values, as defined in Appendix 1
library("parallel")
cl <- makeCluster(getOption ("cl.cores", 2))  	#Can increase the number of clusters

cnt <- 0; 

 
  
  cnt <- cnt + 1
  
  # Generate  x_ij
  lxN <- mulx + rnorm(N)*sqrt(sig2lx)
  lXN <- cbind(rep(1,N), lxN)
  
  ## Generate data from model and select sample
  popllmc <- genpopmixedll.comb.ebdistfuns(D, CNis,edistfun ,e.parms ,bdistfun, sig2lu, lxN, GN, beta0, beta1)
  smc <- popllmc[[4]]; lyNmc <- popllmc[[1]]; yNmc <- popllmc[[2]]; ybarNis <- popllmc[[3]]
  q.pop <- as.vector(popllmc[[5]])	
  lys <- lyNmc[smc]
  lXs <- lXN[smc,]
  u.pop <- popllmc[[6]]
  u.pops <- rbind(u.pops, u.pop)
  dat.temp <- data.frame(Y = lys, X = lXs[,2], area = areafac.pop[smc])

  save.image("TestDataSetNoTrans14March2019.Rdata")

