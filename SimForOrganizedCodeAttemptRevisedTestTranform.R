rm(list=ls(all=TRUE))

##### Set the working directory to a folder that contains the files posted to the Github repository
setwd("G:/Researchers-Investigators/Berg/BaseCodeQRSAE/SortFixBootstrapFromGit9-19-2018/SmallAreaMixedEffectsQuantiles-SortFixBootstrap")

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
Bmax <- 100


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

  ##### ###### ####### ######### ######## REML and EBP ########  ######## ######## ######## ######## ######## 

  ########### Obtain REML estimates of model parameters  
  grmc <- getremlests(smc, D, lyNmc, lXN, GN)
  sig2uhatmc <- grmc[[1]]
  sig2ehatmc <- grmc[[2]]
  betahatmc <- grmc[[3]]
  vhatbetahatmc <- grmc[[4]]

  ##########  Store REML estimates of model parameters
  est.remls <- rbind(est.remls, c(betahatmc, sig2uhatmc, sig2ehatmc))
  
  lbarsi <- t(Gs)%*%lyNmc[smc]/nis
  dbarsi <- t(Gs)%*%lXN[smc,]/nis
  gammais <- sig2uhatmc/(sig2uhatmc + sig2ehatmc/nis)
  mu.dev <- as.vector(gammais*(lbarsi - dbarsi%*%betahatmc))
  mean.cond <- lXN%*%betahatmc + GN%*%(mu.dev) 
  var.cond <- GN%*%(gammais*sig2ehatmc) + sig2ehatmc

  t.min <- min(all.q)
  t.max <- max(all.q)	
  t.seq <- seq(t.min, t.max, length = 100)
  
  ##########  Estimate CDF for conditional normal distribution (not used) 
  cdf.norm.all <- sapply(t.seq, cnp, mean.cond, sqrt(var.cond))
  cdf.smc <- sapply(t.seq, function(t){ ifelse(lys <= t, 1, 0)})
  cdf.norm.all[smc,] <- cdf.smc

  ########## Implement EBP, as described in Remark 2 of Molina and Rao (2010): 
  ########## Simulate a population from the conditional normal distribution (step b of Remark 2)
  r.c.n <- replicate(100, rnorm(length(mean.cond), mean = mean.cond, sd = sqrt(var.cond)))
  ########## Augment the nonsampled units with the vector of sampled elements (step c of Remark 2)
  r.c.n[smc,] <- lys
  ########## Back-tranform the generated population
  r.c.n <- apply(r.c.n, 2, hinvtransG1, lam)

  ########## Compute the quantile for each generated population
  EB.norm <- sapply(1:D, function(d){apply(apply(r.c.n[areafac.pop == d,], 2, quantile, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))[c(1,2,3,4,5),], 1, mean)}) 
  ########## Crude estimate of the variance of the leading term in the EBP
  VEB.norm <- sapply(1:D, function(d){apply(apply(r.c.n[areafac.pop == d,], 2, quantile, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))[c(1,2,3,4,5),]^2, 1, mean)}) - EB.norm^2

  ########## Store EBP estimates of quantiles 
  EB.norm.25s <- rbind(EB.norm.25s, EB.norm[2,])
  EB.norm.50s <- rbind(EB.norm.50s, EB.norm[3,])
  EB.norm.75s <- rbind(EB.norm.75s, EB.norm[4,])
  EB.norm.10s <- rbind(EB.norm.10s, EB.norm[1,])
  EB.norm.90s <- rbind(EB.norm.90s, EB.norm[5,])
  
  ######### Alternative estimator of finite population variance for previous comparison purposes (not used)
  areavar.norm <- sapply(1:D, function(d){ var(as.vector(r.c.n[areafac.pop == d,]))})
  
  ######### Store crude estimates of variances of leading term
  VEB.norm.25s <- rbind(VEB.norm.25s, VEB.norm[2,])
  VEB.norm.50s <- rbind(VEB.norm.50s, VEB.norm[3,])
  VEB.norm.75s <- rbind(VEB.norm.75s, VEB.norm[4,])

  ######### EBP estimate of the finite population mean 
  norm.mean1 <- mean.cond
  norm.mean1[smc] <- lys
  ybar.norms <- rbind(ybar.norms, as.vector(t(GN)%*%norm.mean1)/Nis)

  ############################################     LIGPD Predictors    #################################################################
  
  ########  Compute the LIGPD predictors
  source("estconstrainedtrans_fast_boot.R")
  
  ######## Store the estimates 
  w25JCs <- rbind(w25JCs, areapred[[1]])
  w50JCs <- rbind(w50JCs, areapred[[2]])
  w75JCs <- rbind(w75JCs, areapred[[3]])
  w90JCs <- rbind(w90JCs, areapred[[4]])
  w10JCs <- rbind(w10JCs, areapred[[5]])

  meanJs <- rbind(meanJs, areapred[[6]])
  varJs <- rbind(varJs, areapred[[7]])
  
  mubJs <- rbind(mubJs, areapred[[8]])
  vebJs <- rbind(vebJs, areapred[[9]])

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
    bootsampyxq <- bootgensimple2trans(sig2bhatupdate, GN, XBhatupdate, CNis, nis, D, areafac.pop, rholxilrhouxiu, lamcur)
    q25popbs <- rbind(q25popbs, bootsampyxq[[3]]); q50popbs <- rbind(q50popbs, bootsampyxq[[4]]); q75popbs <- rbind(q75popbs, bootsampyxq[[5]])#
    ########  Construct bootstrap estimates and store bootstrap estimates of the population quantiles
    bootestyxq <- par.bootestsimple2FixSortTransRev(data.frame(Y = bootsampyxq[[1]], X = lxN[bootsampyxq[[2]]]), X = lxN[bootsampyxq[[2]]], areafac.pop[bootsampyxq[[2]]], D, GN[bootsampyxq[[2]],], tauvec,  1500, lxN, b.dist, areafac.pop, bootsampyxq[[2]], use.cl = FALSE, trunc = FALSE, lam = lamcur)
    qhat25bs <- rbind(qhat25bs, bootestyxq[[1]])
    qhat50bs <- rbind(qhat50bs, bootestyxq[[2]])
    qhat75bs <- rbind(qhat75bs, bootestyxq[[3]])
    if(b == Bmax){break}
  }
   MSEhat25b <-  apply((q25popbs - qhat25bs)^2, 2, mean)
   MSEhat50b <-  apply((q50popbs - qhat50bs)^2, 2, mean)
   MSEhat75b <-  apply((q75popbs - qhat75bs)^2, 2, mean)

  mhb25s <- rbind(mhb25s, MSEhat25b)
  mhb50s  <- rbind(mhb50s, MSEhat50b)
  mhb75s  <- rbind(mhb75s, MSEhat75b)
   
  }

  ### Save the Image Every 10 Iterations

  if(cnt%%10 == 0){ save.image("SNTransBoot3Normal.Rdata") }	
  
  print(paste(cnt))
  
  if(cnt == 200){break}
  
}




time.end.all <- Sys.time()



 



