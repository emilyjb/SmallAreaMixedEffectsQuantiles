bootestsimple1FixSort <- function(dat.tempb, X, areafac.samp, D, Gs, tauvec,  Rb, lxN, b.dist, areafac.pop, smc, use.cl, trunc =FALSE){
  
  initpars <- intparfun(dat.tempb, X, areafac.samp,D, Gs,tauvec,  lxN, use.cl)
  
  rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.tempb, tauvec)
  
  XBbetasiguupdate <- updatebetasig2bfunSortFix(tauvec,dat.tempb, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl, Rb )
  
  sig2bhatupdate <- XBbetasiguupdate[[3]]
  betahat <- XBbetasiguupdate[[2]]
  XBhatupdate <- XBbetasiguupdate[[1]]
  
  XBbetasiguupdate2 <- updatebetasig2bfunSortFix(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat , rholxilrhouxiu , areafac.pop, smc, XBhatupdate[smc,], Gs,lxN, use.cl, Rb = 150 )
  
  sig2bhatupdate <- XBbetasiguupdate2[[3]]
  betahat <- XBbetasiguupdate2[[2]]
  XBhatupdate <- XBbetasiguupdate2[[1]]
  
  areapred <- predJW(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl, Rb = 150, trunc )
  
  qhat25b <- areapred[[1]]
  qhat50b <- areapred[[2]]
  qhat75b <- areapred[[3]]
  
  list(qhat25b, qhat50b, qhat75b, sig2bhatupdate, betahat, XBhatupdate, rholxilrhouxiu)
  
}



bootestsimple1FixSortTrans <- function(dat.tempb, X, areafac.samp, D, Gs, tauvec,  Rb, lxN, b.dist, areafac.pop, smc, use.cl, trunc =FALSE, lam){
  
  initpars <- intparfun(dat.tempb, X, areafac.samp,D, Gs,tauvec,  lxN, use.cl)
  
  rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.tempb, tauvec)
  
  XBbetasiguupdate <- updatebetasig2bfunSortFix(tauvec,dat.tempb, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl, Rb )
  
  sig2bhatupdate <- XBbetasiguupdate[[3]]
  betahat <- XBbetasiguupdate[[2]]
  XBhatupdate <- XBbetasiguupdate[[1]]
   
  areapred <- predJWTransFix(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl, Rb = 150, trunc, lam )
  
  qhat25b <- areapred[[1]]
  qhat50b <- areapred[[2]]
  qhat75b <- areapred[[3]]
  
  list(qhat25b, qhat50b, qhat75b, sig2bhatupdate, betahat, XBhatupdate, rholxilrhouxiu)
  
}




bootestsimple2FixSortTrans <- function(dat.tempb, X, areafac.samp, D, Gs, tauvec,  Rb, lxN, b.dist, areafac.pop, smc, use.cl, trunc =FALSE, lam){
  
  initpars <- intparfun(dat.tempb, X, areafac.samp,D, Gs,tauvec,  lxN, use.cl)
  
  rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.tempb, tauvec)
  
  XBbetasiguupdate <- updatebetasig2bfunSortFix(tauvec,dat.tempb, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl, Rb )
  
  sig2bhatupdate <- XBbetasiguupdate[[3]]
  betahat <- XBbetasiguupdate[[2]]
  XBhatupdate <- XBbetasiguupdate[[1]]
   
  areapred <- predJWTransFixB0Exp(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl, Rb = 150, trunc, lam ,0)
  
  qhat25b <- areapred[[1]]
  qhat50b <- areapred[[2]]
  qhat75b <- areapred[[3]]
  
  list(qhat25b, qhat50b, qhat75b, sig2bhatupdate, betahat, XBhatupdate, rholxilrhouxiu)
  
}




bootestsimple2FixSort <- function(dat.tempb, X, areafac.samp, D, Gs, tauvec,  Rb, lxN, b.dist, areafac.pop, smc, use.cl, trunc =FALSE){
  
  initpars <- intparfun(dat.tempb, X, areafac.samp,D, Gs,tauvec,  lxN, use.cl)
  
  rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.tempb, tauvec)
  
  XBbetasiguupdate <- updatebetasig2bfunSortFix(tauvec,dat.tempb, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl, Rb )
  
  sig2bhatupdate <- XBbetasiguupdate[[3]]
  betahat <- XBbetasiguupdate[[2]]
  XBhatupdate <- XBbetasiguupdate[[1]]
   
  areapred <- predJWFix(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl, Rb = 150, trunc )
  
  qhat25b <- areapred[[1]]
  qhat50b <- areapred[[2]]
  qhat75b <- areapred[[3]]
  
  list(qhat25b, qhat50b, qhat75b, sig2bhatupdate, betahat, XBhatupdate)
  
}

