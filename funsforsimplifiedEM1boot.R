bootgensimple1 <- function(sig2bhat.update, GN, XBhat.update1, CNis, nis, D, areafac.pop){

  ahats.b <- rnorm(D, mean = 0, sd = sqrt(sig2bhat.update ))
  all.y.b <- XBhat.update1 + as.vector(GN%*%ahats.b)

  #ran.taus <-  sample(1:length(q.vec), size = nrow(all.y.b), replace = TRUE)
  #mu.boots <- sapply(1:nrow(all.y.b), function(i){ all.y.b[i,ran.taus[i]]})
  
  lyN.b <- apply(all.y.b, 1, sample, size = 1, replace = FALSE)

  #t.min.boot <- min(all.y.b)
  #t.max.boot <- max(all.y.b)
  #t.seq.boot <- seq(t.min.boot, t.max.boot, length = length(q.vec))
  ####  Get the area quantile for the new all.y.b

  #q25b <- sapply(1:D, get.qpop.boot, 0.25, q.vec, all.y.b, t.min.boot, t.max.boot, t.seq.boot)
  #q50b <- sapply(1:D, get.qpop.boot, 0.50, q.vec, all.y.b, t.min.boot, t.max.boot, t.seq.boot)
  #q75b <- sapply(1:D, get.qpop.boot, 0.75, q.vec, all.y.b)

  #sapply(1:D, function(x){sample(1:Nx, size = )})

  smc.b <- sample(1:CNis[1], size=nis[1], replace=FALSE)
  i <- 1
  repeat{
    i <- i + 1
    smc.b <- c(smc.b, sample((CNis[(i-1)] +1):CNis[i], size=nis[i], replace=FALSE))
    if(i == D){break}
  }


  #lyN.b <- apply(all.y.b, 1, sample, size = 1)

  q25boot  <-  tapply(lyN.b, areafac.pop, quantile, prob = 0.25)
  q50boot <-   tapply(lyN.b, areafac.pop, quantile, prob = 0.5)
  q75boot <-  tapply(lyN.b, areafac.pop, quantile, prob = 0.75)

  #q75b <- c(q25b, q5b, q75b)

  lys.b <- lyN.b[smc.b]

  list(lys.b, smc.b, q25boot = q25boot, q50boot = q50boot, q75boot = q75boot)

}

bootestsimple1 <-  function(dat.tempb, X, areafac.samp, D, Gs, tauvec,  Rb, lxN, b.dist, areafac.pop, smc, use.cl, trunc =FALSE){

  initpars <- intparfun(dat.tempb, X, areafac.samp,D, Gs,tauvec,  lxN, use.cl)

  rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.tempb, tauvec)

  XBbetasiguupdate <- updatebetasig2bfun(tauvec,dat.tempb, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl, Rb )

  sig2bhatupdate <- XBbetasiguupdate[[3]]
  betahat <- XBbetasiguupdate[[2]]
  XBhatupdate <- XBbetasiguupdate[[1]]

  XBbetasiguupdate2 <- updatebetasig2bfun(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat , rholxilrhouxiu , areafac.pop, smc, XBhatupdate[smc,], Gs,lxN, use.cl, Rb = 150 )

  sig2bhatupdate <- XBbetasiguupdate2[[3]]
  betahat <- XBbetasiguupdate2[[2]]
  XBhatupdate <- XBbetasiguupdate2[[1]]

  areapred <- predJW(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl, Rb = 150, trunc )

  qhat25b <- areapred[[1]]
  qhat50b <- areapred[[2]]
  qhat75b <- areapred[[3]]
  
  list(qhat25b, qhat50b, qhat75b, sig2bhatupdate, betahat, XBhatupdate)
  
}

  
  