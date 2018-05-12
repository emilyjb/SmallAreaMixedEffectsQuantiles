
rm(list=ls(all=TRUE))
library(nlme)
library("lqmm")
library("sae")
library("survey")
library("quantreg")
library("sn")
library("eva")

source("TwoFunctions.R")
source("gpdfunctionsforwangmethod.R")
source("paroptemfuns.R")
source("EMfunctions2.R")
source("llhoodmaxnorm1.R")
source("genpopmixedllcomb.R")
source("betahatglsfun.R")
source("storoutputsimpleEMsim.R")
source("ebpfuns.R")
source("funsforsimplifiedEM1.R")
source("funsforsimplifiedEM1boot.R")
source("bootgensimple2fun.R")
source("interpkfun.R")
source("updatefix.R")
source("condbfuns1.R")
source("initupdateconstrainedsort.R")
source("estconstrainedfun.R")
library("rmultil")
### Population and sample size configuration:

N <- 20000
B <- 100
#Nis <- c(1200, 4000, 400, 2500, 500, 2000, 800, 7300, 3300, 4000, 3500, 3000, 5000, 2700)
#nis <- 0.1*Nis
nis <- rep(c(rep(5,5),rep(10,5), rep(20,5)),times = 4)
Nis <- round(N*nis/sum(nis),0)
CNis <- cumsum(Nis)
D <- length(nis)
n <- sum(nis)
N <- sum(Nis)
# Sampling fractions
pis  <- nis/Nis
pissamp <- rep(nis/Nis, times = nis)
Nispop <- rep(Nis, Nis)
areafac.pop <- rep(1:D, Nis)
Nissamp <- rep(Nis, nis)

GN <- model.matrix(lm(Nispop~as.factor(areafac.pop)-1))

### Pick a set and set parameters:

beta0 <- -1.5; beta1 <- 0.5
sig2le <- 1 ; sig2lu <- 0.5
mulx <- 0
sig2lx <- 1

M <- 1000; 
cnt <- 0; 

####  Set options:
e.dist <- "Chi"
b.dist <- "Laplace"
doBoot <- TRUE

time.start.all <- Sys.time()

options("contrasts" = c("contr.sum", "contr.poly"))
library("parallel")
### cl <- makeCluster(getOption ("cl.cores", 2))

repeat{
  
  cnt <- cnt + 1
  
  # Generate log(x_ij)
  lxN <- mulx + rnorm(N)*sqrt(sig2lx)
  lXN <- cbind(rep(1,N), lxN)
  
  ## Generate log linear model and select sample
  popllmc <- genpopmixedll.comb.ebdist(D, CNis, sig2le, sig2lu,beta0, beta1, Nis, lxN, GN,   2, e.dist, b.dist, 0.75)
  smc <- popllmc[[4]]; lyNmc <- popllmc[[1]]; yNmc <- popllmc[[2]]; ybarNis <- popllmc[[3]]
  q.pop <- as.vector(popllmc[[5]])	
  lys <- lyNmc[smc]
  lXs <- lXN[smc,]
  u.pop <- popllmc[[6]]
  u.pops <- rbind(u.pops, u.pop)
  dat.temp <- data.frame(Y = lys, X = lXs[,2], area = areafac.pop[smc])
  plot(dat.temp[,c(2,1)])
  ##### ###### ####### ######### ########  ALD Procedure  ########  ######## ######## ######## ######## ######## 
  qVecAld  <- seq(0.01, 0.99, by = 0.01)
  lq.fit <- lqmm(Y~X, data = dat.temp, random = ~1, tau = qVecAld, group = area)
  alpha.hats <- ranef.lqmm(lq.fit)
  ahat.mat <- matrix(unlist( ranef.lqmm(lq.fit)), nrow = length(unique(areafac.pop)), byrow = FALSE)
  beta.hats <-  coef.lqmm(lq.fit)
  bhat.mat <- matrix(unlist( coef.lqmm(lq.fit)), nrow = 2, byrow = FALSE)
  scale.lqs <- c(scale.lqs, VarCorr(lq.fit)[[which(qVecAld == 0.5)]])
  
  b0als <- rbind(b0als, bhat.mat[1,])
  b1als <- rbind(b1als, bhat.mat[2,])
  aials <- rbind(aials, apply(ahat.mat, 1, mean))

  all.q <- lXN%*%bhat.mat + GN%*%ahat.mat

  #q.MC <- t(apply(all.q, 1, sample, size = 1000, replace = TRUE))
  
  #w25s <- rbind(w25s, sapply(1:D, function(i){ quantile(as.vector( q.MC[areafac.pop == i,]), prob = 0.25)}))
  #w50s <- rbind(w50s, sapply(1:D, function(i){ quantile(as.vector( q.MC[areafac.pop == i,]), prob = 0.5)}))
  #w75s <- rbind(w75s, sapply(1:D, function(i){ quantile(as.vector( q.MC[areafac.pop == i,]), prob = 0.75)}))
  
  q25bs <- rbind(q25bs, sapply(1:D, function(i){ quantile(as.vector( all.q[areafac.pop == i,]), prob = 0.25)}))
  q5bs <- rbind(q5bs, sapply(1:D, function(i){ quantile(as.vector( all.q[areafac.pop == i,]), prob = 0.5)}))
  q75bs <- rbind(q75bs, sapply(1:D, function(i){ quantile(as.vector( all.q[areafac.pop == i,]), prob = 0.75)}))
  q10bs <- rbind(q10bs, sapply(1:D, function(i){ quantile(as.vector( all.q[areafac.pop == i,]), prob = 0.1)}))
  q90bs <- rbind(q90bs, sapply(1:D, function(i){ quantile(as.vector( all.q[areafac.pop == i,]), prob = 0.9)}))
  
  meanalss <- rbind(meanalss,  sapply(1:D, function(i){ mean(as.vector( all.q[areafac.pop == i,]))}))
  varalss <- rbind(varalss,  sapply(1:D, function(i){ var(as.vector( all.q[areafac.pop == i,]))}))
  
  ##### ###### ####### ######### ########  Sample Quantiles and Means ########  ######## ######## ######## ######## ######## 
  
  qhats.samp10 <- rbind(qhats.samp10,  tapply(lys, areafac.pop[smc], quantile, prob = 0.10))
  qhats.samp25 <- rbind(qhats.samp25,  tapply(lys, areafac.pop[smc], quantile, prob = 0.25))
  qhats.samp5 <- rbind(qhats.samp5,  tapply(lys, areafac.pop[smc], quantile, prob = 0.5))
  qhats.samp75 <- rbind(qhats.samp75,  tapply(lys, areafac.pop[smc], quantile, prob = 0.75))
  qhats.samp90 <- rbind(qhats.samp75,  tapply(lys, areafac.pop[smc], quantile, prob = 0.90))
  
  Gs <- GN[smc,] 
  ybarnis <- t(Gs)%*%lyNmc[smc]/nis
  lybarniss <- rbind(lybarniss, ybarnis)
  
  ##### ###### ####### ######### ######## Population Quantiles and Means  ########  ######## ######## ######## ######## ######## 
  
  qhats.pop10 <- rbind(qhats.pop10, tapply(lyNmc , areafac.pop, quantile, prob = 0.1))
  qhats.pop25 <- rbind(qhats.pop25, tapply(lyNmc , areafac.pop, quantile, prob = 0.25))
  qhats.pop5 <- rbind(qhats.pop5, tapply(lyNmc , areafac.pop, quantile, prob = 0.5))
  qhats.pop75 <- rbind(qhats.pop75, tapply(lyNmc , areafac.pop, quantile, prob = 0.75))
  qhats.pop90 <- rbind(qhats.pop90, tapply(lyNmc , areafac.pop, quantile, prob = 0.90))

  lbarNis <- as.vector(t(GN)%*%lyNmc)/Nis
  lbarNiss <- rbind(lbarNiss, lbarNis)
  
  ##### ###### ####### ######### ######## REML and EBP ########  ######## ######## ######## ######## ######## 
  
  grmc <- getremlests(smc, D, lyNmc, lXN, GN)
  sig2uhatmc <- grmc[[1]]
  sig2ehatmc <- grmc[[2]]
  betahatmc <- grmc[[3]]
  vhatbetahatmc <- grmc[[4]]
  
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
  
  cdf.norm.all <- sapply(t.seq, cnp, mean.cond, sqrt(var.cond))
  cdf.smc <- sapply(t.seq, function(t){ ifelse(lys <= t, 1, 0)})
  cdf.norm.all[smc,] <- cdf.smc

  r.c.n <- replicate(100, rnorm(length(mean.cond), mean = mean.cond, sd = sqrt(var.cond)))
  r.c.n[smc,] <- lys
  EB.norm <- sapply(1:D, function(d){apply(apply(r.c.n[areafac.pop == d,], 2, quantile, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))[c(1,2,3,4,5),], 1, mean)})
  
  VEB.norm <- sapply(1:D, function(d){apply(apply(r.c.n[areafac.pop == d,], 2, quantile, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))[c(1,2,3,4,5),]^2, 1, mean)}) - EB.norm^2
  
  EB.norm.25s <- rbind(EB.norm.25s, EB.norm[2,])
  EB.norm.50s <- rbind(EB.norm.50s, EB.norm[3,])
  EB.norm.75s <- rbind(EB.norm.75s, EB.norm[4,])
  EB.norm.10s <- rbind(EB.norm.10s, EB.norm[1,])
  EB.norm.90s <- rbind(EB.norm.90s, EB.norm[5,])
  
  areavar.norm <- sapply(1:D, function(d){ var(as.vector(r.c.n[areafac.pop == d,]))})
  
  VEB.norm.25s <- rbind(VEB.norm.25s, VEB.norm[1,])
  VEB.norm.50s <- rbind(VEB.norm.50s, VEB.norm[2,])
  VEB.norm.75s <- rbind(VEB.norm.75s, VEB.norm[3,])
  
  norm.mean1 <- mean.cond
  norm.mean1[smc] <- lys
  ybar.norms <- rbind(ybar.norms, as.vector(t(GN)%*%norm.mean1)/Nis)

  ############################################     JW Predictors    #################################################################
  
  source("estconstrainedfix.R")
 
 # w25Js <- rbind(w25Js, areapred[[1]])
 # w50Js <- rbind(w50Js, areapred[[2]])
 # w75Js <- rbind(w75Js, areapred[[3]])
 # w90Js <- rbind(w90Js, areapred[[4]])
 # w10Js <- rbind(w10Js, areapred[[5]])
  
  meanJs <- rbind(meanJs, areapred[[6]])
  varJs <- rbind(varJs, areapred[[7]])
  
  mubJs <- rbind(mubJs, areapred[[8]])
  vebJs <- rbind(vebJs, areapred[[9]])
  
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
    bootsampyxq <- bootgensimple3(sig2bhatupdate, GN, XBhatupdate, CNis, nis, D, areafac.pop, rholxilrhouxiu, smc, b.dist)
    q25popbs <- rbind(q25popbs, bootsampyxq[[3]]); q50popbs <- rbind(q50popbs, bootsampyxq[[4]]); q75popbs <- rbind(q75popbs, bootsampyxq[[5]])
    bootestyxq <- bootestsimple1(data.frame(Y = bootsampyxq[[1]], X = lxN[bootsampyxq[[2]]]), X = lxN[bootsampyxq[[2]]], areafac.pop[bootsampyxq[[2]]], D, GN[bootsampyxq[[2]],], tauvec,  1500, lxN, b.dist, areafac.pop, bootsampyxq[[2]], use.cl = FALSE)
    qhat25bs <- rbind(qhat25bs, bootestyxq[[1]])
    qhat50bs <- rbind(qhat50bs, bootestyxq[[2]])
    qhat75bs <- rbind(qhat75bs, bootestyxq[[3]])
    if(b == B){break}
  }
   MSEhat25b <-  apply((q25popbs - qhat25bs)^2, 2, mean)
   MSEhat50b <-  apply((q50popbs - qhat50bs)^2, 2, mean)
   MSEhat75b <-  apply((q75popbs - qhat75bs)^2, 2, mean)

   mhb25s <- rbind(mhb25s, MSEhat25b)
   mhb50s  <- rbind(mhb50s, MSEhat50b)
   mhb75s  <- rbind(mhb75s, MSEhat75b)
   
  }

  if(cnt%%10 == 0){ save.image("ChiOutput/ChiOutputInterpF1LaplaceFix.Rdata") }	
  
  print(paste(cnt))
  
  if(cnt == 200){break}
  
}




time.end.all <- Sys.time()


(mean(apply( mhb25s, 2, mean)) - mean(apply((qhats.pop25  - w25JCs)^2, 2, mean)))/ mean(apply((qhats.pop25  - w25JCs)^2, 2, mean))*100






