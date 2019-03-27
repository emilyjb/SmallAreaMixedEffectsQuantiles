rm(list = ls(all = TRUE))

#####  Load test data set:

load("TestDataSetNoTrans24March2019.Rdata" )


###  Load libraries
library(nlme)
library("lqmm")
library("sae")
library("survey")
library("quantreg")
library("sn")
library("eva")
library("rmutil")

options("contrasts" = c("contr.sum", "contr.poly"))  ### Use sum to zero constraints for initial values, as defined in Appendix 1
library("parallel")
cl <- makeCluster(getOption ("cl.cores", 2))  	#Can increase the number of clusters


#####  Compute ALD, NEB, and LIGPD predictors for test data set:

 qVecAld  <- seq(0.01, 0.99, by = 0.01)

  ###########  Estimate parameters of ALD procedure
  lq.fit <- lqmm(Y~X, data = dat.temp, random = ~1, tau = qVecAld, group = area, nK = 30, type = aldtype)
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
  qhats.samp10 <- rbind(qhats.samp10,  tapply(lys, areafac.pop[smc], quantile, prob = 0.10))
  qhats.samp25 <- rbind(qhats.samp25,  tapply(lys, areafac.pop[smc], quantile, prob = 0.25))
  qhats.samp5 <- rbind(qhats.samp5,  tapply(lys, areafac.pop[smc], quantile, prob = 0.5))
  qhats.samp75 <- rbind(qhats.samp75,  tapply(lys, areafac.pop[smc], quantile, prob = 0.75))
  qhats.samp90 <- rbind(qhats.samp90,  tapply(lys, areafac.pop[smc], quantile, prob = 0.90))

  ########### Calculate and store sample means  
  Gs <- GN[smc,] 
  ybarnis <- t(Gs)%*%lyNmc[smc]/nis
  lybarniss <- rbind(lybarniss, ybarnis)
    
  ##### ###### ####### ######### ######## Population Quantiles and Means  ########  ######## ######## ######## ######## ######## 

  ########### Calculate and store population quantiles
  qhats.pop10 <- rbind(qhats.pop10, tapply(lyNmc , areafac.pop, quantile, prob = 0.1))
  qhats.pop25 <- rbind(qhats.pop25, tapply(lyNmc , areafac.pop, quantile, prob = 0.25))
  qhats.pop5 <- rbind(qhats.pop5, tapply(lyNmc , areafac.pop, quantile, prob = 0.5))
  qhats.pop75 <- rbind(qhats.pop75, tapply(lyNmc , areafac.pop, quantile, prob = 0.75))
  qhats.pop90 <- rbind(qhats.pop90, tapply(lyNmc , areafac.pop, quantile, prob = 0.90))

  ########### Calculate and store population means
  lbarNis <- as.vector(t(GN)%*%lyNmc)/Nis
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
  var.cond <- GN%*%(gammais*sig2ehatmc/nis) + sig2ehatmc

  mudevs <- rbind(mudevs, mu.dev)


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
  source("estconstrainedfix_fast_boot.R")
  
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



  #####  Summarize distributions of estimates of quantiles for comparison with Github documentation:

  ################## ALD estimates
  summary(q10bs[1,])
  summary(q25bs[1,])
  summary(q5bs[1,])
  summary(q75bs[1,])
  summary(q90bs[1,])

  ################## Normal EB estimates
  summary(EB.norm.10s[1,])
  summary(EB.norm.25s[1,])
  summary(EB.norm.50s[1,])
  summary(EB.norm.75s[1,])
  summary(EB.norm.90s[1,])

  ################## LIGPD estimates
  summary(w10JCs[1,])
  summary(w25JCs[1,])
  summary(w50JCs[1,])
  summary(w75JCs[1,])
  summary(w90JCs[1,])

 ####  Summarize in data set:
 sumests <- rbind( summary(q10bs[1,]),
  summary(q25bs[1,]),
  summary(q5bs[1,]),
  summary(q75bs[1,]),
  summary(q90bs[1,]),

  summary(EB.norm.10s[1,]),
  summary(EB.norm.25s[1,]),
  summary(EB.norm.50s[1,]),
  summary(EB.norm.75s[1,]),
  summary(EB.norm.90s[1,]),

  summary(w10JCs[1,]),
  summary(w25JCs[1,]),
  summary(w50JCs[1,]),
  summary(w75JCs[1,]),
  summary(w90JCs[1,])  
)

 data.frame( Method = rep(c("ALD", "NEB", "LIGPD"), each = 5),QuantileLevel = rep(c(10, 25, 50, 75, 90), times = 3), sumests)
 
 dfout <- data.frame( Method = rep(c("ALD", "NEB", "LIGPD"), each = 5),QuantileLevel = rep(c(10, 25, 50, 75, 90), times = 3), sumests)






