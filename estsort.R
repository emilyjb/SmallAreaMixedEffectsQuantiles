
tauvec <- seq(0.01, 0.99, by = 0.01)
areafac.samp <- areafac.pop[smc]
Gs <- GN[smc,]
X <- lxN[smc]

initpars <-  intparfun(dat.temp, X, areafac.samp,D, Gs,tauvec, lxN, use.cl = FALSE)
 
rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.temp, tauvec)

XBbetasiguupdate <- par.updatebetasig2bfunSortFix (tauvec,dat.temp, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl= FALSE, Rb = 1500 )


sig2bhatupdate <- XBbetasiguupdate[[3]]
betahat <- XBbetasiguupdate[[2]]
XBhatupdate <- XBbetasiguupdate[[1]]

rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.temp, tauvec)


areapred  <- par.predJWFix(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl= FALSE, Rb = 1500 )

  w25JSs <- rbind(w25JSs, areapred[[1]])
  w50JSs <- rbind(w50JSs, areapred[[2]])
  w75JSs <- rbind(w75JSs, areapred[[3]])
  w90JSs <- rbind(w90JSs, areapred[[4]])
  w10JSs <- rbind(w10JSs, areapred[[5]])



XBbetasiguupdate2 <- par.updatebetasig2bfunSortFix(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat , rholxilrhouxiu , areafac.pop, smc, XBhatupdate[smc,], Gs,lxN, use.cl= FALSE, Rb = 1500 )

sig2bhatupdate <- XBbetasiguupdate2[[3]]
betahat <- XBbetasiguupdate2[[2]]
XBhatupdate <- XBbetasiguupdate2[[1]]

bj10Ss <- rbind(bj10Ss, betahat[1,])#
bj11Ss <- rbind(bj11Ss, betahat[2,])

rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.temp, tauvec)


areapred <- par.predJWFix(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl= FALSE, Rb = 1500 )


  w25JS2s <- rbind(w25JS2s, areapred[[1]])
  w50JS2s <- rbind(w50JS2s, areapred[[2]])
  w75JS2s <- rbind(w75JS2s, areapred[[3]])
  w90JS2s <- rbind(w90JS2s, areapred[[4]])
  w10JS2s <- rbind(w10JS2s, areapred[[5]])



