rm(list=ls(all=TRUE))

##### Set the working directory to a folder that contains the files posted to the Github repository
setwd("C:/Users/Emily/Documents/GitHub/Simulations/SmallAreaMixedEffectsQuantiles-CompareExpBoToBoExpPluginOpt")

### Load libraries
library(nlme)
library("lqmm")
library("sae")
library("survey")
library("quantreg")
library("sn")
library("eva")


## Source R files
source("storoutputsimpleEMsim.R")		# Store output
source("genpopmixedllcomb.R")			# Generate data
source("transq1.R")				# For the transformation
source("NEBfuns.R")				# For EBP predictors
source("LIGPDfuns.R")			# For LIGPD predictors 
source("LIGPDfuns_trans.R")			# For LIGPD predictors - Transformation
source("Bootfuns.R")				# For Bootstrap MSE computation


### Population and sample size configuration:
N <- 20000
nis <- rep(c(rep(5,5),rep(10,5), rep(20,5)),times = 4)
Nis <- round(N*nis/sum(nis),0)
CNis <- cumsum(Nis)
D <- length(nis)
n <- sum(nis)
N <- sum(Nis)

### Bootstrap sample size (denoted by T in manuscript)
B  <- 100


# Sampling fractions
pis  <- nis/Nis
pissamp <- rep(nis/Nis, times = nis)
Nispop <- rep(Nis, Nis)
areafac.pop <- rep(1:D, Nis)
Nissamp <- rep(Nis, nis)

GN <- model.matrix(lm(Nispop~as.factor(areafac.pop)-1))

### Set parameters:
beta0 <- -1.5; beta1 <- 0.5
sig2le <- 1; sig2lu <- 0.5
mulx <- 0
sig2lx <- 1

cnt <- 0
lam <- 0

####  Set error and random effect distribution options:
####### Error distribution
e.dist <- "SN"  ### Other options are "Chi" or "T"
####### Area random effect distribution
b.dist <- "Normal"  ###  Other option is "Laplace"
doBoot <- TRUE  #### Change to FALSE to skip bootstrap

time.start.all <- Sys.time()

options("contrasts" = c("contr.sum", "contr.poly"))
library("parallel")
cl <- makeCluster(getOption ("cl.cores", 2))

repeat{
  cnt <- cnt + 1
  
  # Generate  x_ij
  lxN <- mulx + rnorm(N)*sqrt(sig2lx)
  lXN <- cbind(rep(1,N), lxN)
  
  ## Generate data from model and select sample
  popllmc <- genpopmixedll.comb.ebdist(D, CNis, sig2le, sig2lu,beta0, beta1, Nis, lxN, GN,   5, e.dist, b.dist, 0.75)
  smc <- popllmc[[4]]; lyNmc <- popllmc[[1]]; yNmc <- popllmc[[2]]; ybarNis <- popllmc[[3]]
  q.pop <- as.vector(popllmc[[5]])	
  lys <- lyNmc[smc]
  lXs <- lXN[smc,]
  u.pop <- popllmc[[6]]
  u.pops <- rbind(u.pops, u.pop)
  dat.temp <- data.frame(Y = lys, X = lXs[,2], area = areafac.pop[smc])
  
  ##### ###### ####### ######### ########  ALD Procedure  ########  ######## ######## ######## ######## ######## 
  qVecAld  <- seq(0.01, 0.99, by = 0.01)

  ###########  Estimate parameters of ALD procedure
  lq.fit <- lqmm(Y~X, data = dat.temp, random = ~1, tau = qVecAld, group = area)
  ###########  ALD predictors of random effects
  ahat.mat <- matrix(unlist( ranef.lqmm(lq.fit)), nrow = length(unique(areafac.pop)), byrow = FALSE)
  ###########  ALD estimators of coefficients
  bhat.mat <- matrix(unlist( coef.lqmm(lq.fit)), nrow = 2, byrow = FALSE)
  ###########  ALD estimators of between-area variance components
  scale.lqs <- c(scale.lqs, VarCorr(lq.fit)[[which(qVecAld == 0.5)]])

  ###########  Store ALD estimators of regression coefficients and between-area variance
  ######################  ALD estimators of intercept  
  b0als <- rbind(b0als, bhat.mat[1,])
  ######################  ALD estimators of slope
  b1als <- rbind(b1als, bhat.mat[2,])
  ######################  ALD predictors of random effects
  aials <- rbind(aials, apply(ahat.mat, 1, mean))

  ###########  Construct ALD approximation for the population
  all.q <- lXN%*%bhat.mat + GN%*%ahat.mat
  all.q <- apply(all.q, 2, hinvtransG1, lam)

  ###########  Calculate and store ALD predictors of quantiles 
  q25bs <- rbind(q25bs, sapply(1:D, function(i){ quantile(as.vector( all.q[areafac.pop == i,]), prob = 0.25)}))
  q5bs <- rbind(q5bs, sapply(1:D, function(i){ quantile(as.vector( all.q[areafac.pop == i,]), prob = 0.5)}))
  q75bs <- rbind(q75bs, sapply(1:D, function(i){ quantile(as.vector( all.q[areafac.pop == i,]), prob = 0.75)}))
  q10bs <- rbind(q10bs, sapply(1:D, function(i){ quantile(as.vector( all.q[areafac.pop == i,]), prob = 0.1)}))
  q90bs <- rbind(q90bs, sapply(1:D, function(i){ quantile(as.vector( all.q[areafac.pop == i,]), prob = 0.9)}))
 
  ###########  Calculate and store ALD predictors of finite population mean and variance (not in paper) 
  meanalss <- rbind(meanalss,  sapply(1:D, function(i){ mean(as.vector( all.q[areafac.pop == i,]))}))
  varalss <- rbind(varalss,  sapply(1:D, function(i){ var(as.vector( all.q[areafac.pop == i,]))}))
  
  ##### ###### ####### ######### ########  Sample Quantiles and Means ########  ######## ######## ######## ######## ######## 
  
  ########### Calculate and store sample quantiles
  qhats.samp10 <- rbind(qhats.samp10,  tapply(exp(lys), areafac.pop[smc], quantile, prob = 0.10))
  qhats.samp25 <- rbind(qhats.samp25,  tapply(exp(lys), areafac.pop[smc], quantile, prob = 0.25))
  qhats.samp5 <- rbind(qhats.samp5,  tapply(exp(lys), areafac.pop[smc], quantile, prob = 0.5))
  qhats.samp75 <- rbind(qhats.samp75,  tapply(exp(lys), areafac.pop[smc], quantile, prob = 0.75))
  qhats.samp90 <- rbind(qhats.samp90,  tapply(exp(lys), areafac.pop[smc], quantile, prob = 0.90))
 
  ########### Calculate and store sample means   
  Gs <- GN[smc,] 
  ybarnis <- t(Gs)%*%lyNmc[smc]/nis
  lybarniss <- rbind(lybarniss, ybarnis)

  ########################   Note: Sample Quantile and mean invariant to monotone transformation, and back-transformation takes place in output file.
  
  ##### ###### ####### ######### ######## Population Quantiles and Means  ########  ######## ######## ######## ######## ######## 
  
  ########### Calculate and store population quantiles
  qhats.pop10 <- rbind(qhats.pop10, tapply(hinvtransG1(lyNmc, lam) , areafac.pop, quantile, prob = 0.1))
  qhats.pop25 <- rbind(qhats.pop25, tapply(hinvtransG1(lyNmc, lam), areafac.pop, quantile, prob = 0.25))
  qhats.pop5 <- rbind(qhats.pop5, tapply(hinvtransG1(lyNmc, lam), areafac.pop, quantile, prob = 0.5))
  qhats.pop75 <- rbind(qhats.pop75, tapply(hinvtransG1(lyNmc, lam) , areafac.pop, quantile, prob = 0.75))
  qhats.pop90 <- rbind(qhats.pop90, tapply(hinvtransG1(lyNmc, lam), areafac.pop, quantile, prob = 0.90))

  ########### Calculate and store population means
  lbarNis <- as.vector(t(GN)%*%hinvtransG1(lyNmc, lam))/Nis
  lbarNiss <- rbind(lbarNiss, lbarNis)

 
  ############################################     LIGPD Predictors    #################################################################
  
  source("estconstrainedfix_fast_boot.R")

  ######## Store the estimates 
  w25JRs <- rbind(w25JRs, exp(areapred[[1]]))
  w50JRs <- rbind(w50JRs, exp(areapred[[2]]))
  w75JRs <- rbind(w75JRs, exp(areapred[[3]]))
  w90JRs <- rbind(w90JRs, exp(areapred[[4]]))
  w10JRs <- rbind(w10JRs, exp(areapred[[5]]))

  ######## Store the raw estimates 
  w25JCs <- rbind(w25JCs, areapred[[1]])
  w50JCs <- rbind(w50JCs, areapred[[2]])
  w75JCs <- rbind(w75JCs, areapred[[3]])
  w90JCs <- rbind(w90JCs, areapred[[4]])
  w10JCs <- rbind(w10JCs, areapred[[5]])

 # meanJs <- rbind(meanJs, areapred[[6]])
 # varJs <- rbind(varJs, areapred[[7]])
  
 # mubJs <- rbind(mubJs, areapred[[8]])
 # vebJs <- rbind(vebJs, areapred[[9]])


##############################################  Implement bootstrap if doBoot is TRUE ############################################## 
if(doBoot){
  
  q25popbs <- c()
  q50popbs <- c()
  q75popbs <- c()
  qhat25bs <- c()
  qhat50bs <- c()
  qhat75bs <- c()


  b <- 0
  repeat{
    b <- b + 1
    ######## Generate a bootstrap sample and store bootstrap versions of the population quantiles
    bootsampyxq <- bootgensimple3(sig2bhatupdate, GN, XBhatupdate, CNis, nis, D, areafac.pop, rholxilrhouxiu, smc, b.dist)
    q25popbs <- rbind(q25popbs, bootsampyxq[[3]]); q50popbs <- rbind(q50popbs, bootsampyxq[[4]]); q75popbs <- rbind(q75popbs, bootsampyxq[[5]])
    ########  Construct bootstrap estimates and store bootstrap estimates of the population quantiles
    bootestyxq <- par.bootestsimple1Rev(data.frame(Y = bootsampyxq[[1]], X = lxN[bootsampyxq[[2]]]), X = lxN[bootsampyxq[[2]]], areafac.pop[bootsampyxq[[2]]], D, GN[bootsampyxq[[2]],], tauvec,  1500, lxN, b.dist, areafac.pop, bootsampyxq[[2]], use.cl = FALSE)
    qhat25bs <- rbind(qhat25bs, bootestyxq[[1]])
    qhat50bs <- rbind(qhat50bs, bootestyxq[[2]])
    qhat75bs <- rbind(qhat75bs, bootestyxq[[3]])
    if(b == B){break}
  }

   MSEhat25b <-  apply((q25popbs - qhat25bs)^2, 2, mean)
   MSEhat50b <-  apply((q50popbs - qhat50bs)^2, 2, mean)
   MSEhat75b <-  apply((q75popbs - qhat75bs)^2, 2, mean)

   MSEhat25bExp <-  apply((exp(q25popbs) -exp(qhat25bs))^2, 2, mean)
   MSEhat50bExp <-  apply((exp(q50popbs) - exp(qhat50bs))^2, 2, mean)
   MSEhat75bExp <-  apply((exp(q75popbs) - exp(qhat75bs))^2, 2, mean)

   mhb25s <- rbind(mhb25s, MSEhat25b)
   mhb50s  <- rbind(mhb50s, MSEhat50b)
   mhb75s  <- rbind(mhb75s, MSEhat75b)

   mhb25Exps <- rbind(mhb25Exps, MSEhat25bExp)
   mhb50Exps  <- rbind(mhb50Exps, MSEhat50bExp)
   mhb75Exps  <- rbind(mhb75Exps, MSEhat75bExp)



  lci25Exps <- exp(areapred[[1]] - 1.96*sqrt(mhb25s))
  uci25Exps <- exp(areapred[[1]] + 1.96*sqrt(mhb25s))

  lci50Exps <- exp(areapred[[2]] - 1.96*sqrt(mhb50s))
  uci50Exps <- exp(areapred[[2]] + 1.96*sqrt(mhb50s))

  lci75Exps <- exp(areapred[[3]] - 1.96*sqrt(mhb75s))
  uci75Exps <- exp(areapred[[3]] + 1.96*sqrt(mhb75s))


   
  }
  

  ### Save the Image Every 10 Iterations

  if(cnt%%10 == 0){ save.image("EvaluatePluginBoot/PluginBootExp1.Rdata") }	
  
  print(paste(cnt))
  
  if(cnt == 200){break}
  
}




time.end.all <- Sys.time()



 



