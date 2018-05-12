
tauvec <- seq(0.01, 0.99, by = 0.01)
areafac.samp <- areafac.pop[smc]
Gs <- GN[smc,]
X <- lxN[smc]

initpars <- intparfun(dat.temp, X, areafac.samp,D, Gs,tauvec, lxN, use.cl = FALSE)

bj00Cs <- rbind(bj00Cs, initpars$beta[1,])
bj01Cs <- rbind(bj01Cs, initpars$beta[2,])

rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.temp, tauvec)

XBbetasiguupdate <- updatebetasig2bfunConstrFix(tauvec,dat.temp, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl= FALSE, Rb = 1500 )

sig2bhatupdate <- XBbetasiguupdate[[3]]
betahat <- XBbetasiguupdate[[2]]
XBhatupdate <- XBbetasiguupdate[[1]]
 
bj10Cs <- rbind(bj10Cs, betahat[1,])
bj11Cs <- rbind(bj11Cs, betahat[2,])

areapred <- predJWFix(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl= FALSE, Rb = 1500 )


  w25JCs <- rbind(w25JCs, areapred[[1]])
  w50JCs <- rbind(w50JCs, areapred[[2]])
  w75JCs <- rbind(w75JCs, areapred[[3]])
  w90JCs <- rbind(w90JCs, areapred[[4]])
  w10JCs <- rbind(w10JCs, areapred[[5]])


#ystarMC <-replicate(500, areapredraninterp(areapred$all.q.J, rholxilrhouxiu, tauvec))
#ystarMC[smc,] <- lys
#w10JRs <- rbind(w10JRs, sapply(1:D, function(i){ quantile(as.vector( ystarMC[areafac.pop == i,]), prob = 0.1)}))
#w25JRs <- rbind(w25JRs, sapply(1:D, function(i){ quantile(as.vector( ystarMC[areafac.pop == i,]), prob = 0.25)}))
#w50JRs <- rbind(w50JRs, sapply(1:D, function(i){ quantile(as.vector( ystarMC[areafac.pop == i,]), prob = 0.5)}))
#w75JRs <- rbind(w75JRs, sapply(1:D, function(i){ quantile(as.vector( ystarMC[areafac.pop == i,]), prob = 0.75)}))
#w90JRs <- rbind(w90JRs, sapply(1:D, function(i){ quantile(as.vector( ystarMC[areafac.pop == i,]), prob = 0.9)}))
  




 