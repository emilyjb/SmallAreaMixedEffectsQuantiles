
tauvec <- seq(0.01, 0.99, by = 0.01)
areafac.samp <- areafac.pop[smc]
Gs <- GN[smc,]
X <- lxN[smc]

initpars <-  intparfun(dat.temp, X, areafac.samp,D, Gs,tauvec, lxN, use.cl = FALSE)
 
rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.temp, tauvec)

XBbetasiguupdate <- par.updatebetasig2bfunIsoregFix (tauvec,dat.temp, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl= FALSE, Rb = 1500 )


sig2bhatupdate <- XBbetasiguupdate[[3]]
betahat <- XBbetasiguupdate[[2]]
XBhatupdate <- XBbetasiguupdate[[1]]

rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.temp, tauvec)


areapred  <- par.predJWFix(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl= FALSE, Rb = 1500 )

  w25JRs <- rbind(w25JRs, areapred[[1]])
  w50JRs <- rbind(w50JRs, areapred[[2]])
  w75JRs <- rbind(w75JRs, areapred[[3]])
  w90JRs <- rbind(w90JRs, areapred[[4]])
  w10JRs <- rbind(w10JRs, areapred[[5]])


XBbetasiguupdate2 <- par.updatebetasig2bfunIsoregFix(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat , rholxilrhouxiu , areafac.pop, smc, XBhatupdate[smc,], Gs,lxN, use.cl= FALSE, Rb = 1500 )

sig2bhatupdate <- XBbetasiguupdate2[[3]]
betahat <- XBbetasiguupdate2[[2]]
XBhatupdate <- XBbetasiguupdate2[[1]]

bj10Ss <- rbind(bj10Ss, betahat[1,])#
bj11Ss <- rbind(bj11Ss, betahat[2,])

rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.temp, tauvec)


areapred <- par.predJWFix(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl= FALSE, Rb = 1500 )


  w25JR2s <- rbind(w25JR2s, areapred[[1]])
  w50JR2s <- rbind(w50JR2s, areapred[[2]])
  w75JR2s <- rbind(w75JR2s, areapred[[3]])
  w90JR2s <- rbind(w90JR2s, areapred[[4]])
  w10JR2s <- rbind(w10JR2s, areapred[[5]])


