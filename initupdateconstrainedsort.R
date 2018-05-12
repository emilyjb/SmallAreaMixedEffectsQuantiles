intparfunConstr <- function(dat.temp, X, areafac.samp, D, Gs, tauvec, lxN, use.cl = TRUE){
  
  lm1med <- rq(dat.temp$Y~X + as.factor(areafac.samp), tau = 0.5)
  bhats <-  summary.rq(lm1med, se = "ker")$coefficients[,"Value"][-c(1,2)]
  vbhat <-  summary.rq(lm1med, se = "ker")$coefficients[,"Std. Error"][-c(1,2)]^2
  sig2bhat <- betahat.glsfun(vbhat, bhats , 2, matrix(1,nrow = D-1), bhats)[[2]]

  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%c( bhats,  -sum(bhats)))
  coefinit <- estconstrained(dat.temp, as.formula("dev~X"),   tauvec)
  paravec <- c(as.vector(coefinit), c( bhats,  -sum(bhats)))
  XBhat.init <- cbind(1, lxN)%*%coefinit
  list(beta = coefinit, sig2bhat= sig2bhat, bhats = bhats, XBinit = XBhat.init)

}


intparfunSort <- function(dat.temp, X, areafac.samp, D, Gs, tauvec, lxN, use.cl = TRUE){
  
  lm1med <- rq(dat.temp$Y~X + as.factor(areafac.samp), tau = 0.5)
  bhats <-  summary.rq(lm1med, se = "ker")$coefficients[,"Value"][-c(1,2)]
  vbhat <-  summary.rq(lm1med, se = "ker")$coefficients[,"Std. Error"][-c(1,2)]^2
  sig2bhat <- betahat.glsfun(vbhat, bhats , 2, matrix(1,nrow = D-1), bhats)[[2]]

  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%c( bhats,  -sum(bhats)))
  fit.init <- rq(dev~X, data = dat.temp, tau = tauvec)
  paravec <- c(as.vector(fit.init$coef), c( bhats,  -sum(bhats)))
  XBhat.init <- cbind(1, lxN)%*%fit.init$coef
  XBhat.init1 <- t( apply( XBhat.init, 1, sort))
 
  list(beta = fit.init$coef, sig2bhat= sig2bhat, bhats = bhats, XBinit = XBhat.init1)


}


updatebetasig2bfunConstr <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
  seq.points <- qnorm(tauvec, mean = 0, sd = sqrt(sig2bhat))
  lys <- dat.temp$Y
  rhohat.l <- gpdpar[1]; xi.l <- gpdpar[2]; rhohat.u <- gpdpar[3]; xi.u <- gpdpar[4]
  if(use.cl){
    if(bdist == "Laplace"){
      clusterExport(cl, c("seq.points", "lys", "XB.init", "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "compfnum.mb1forintplaplace", "comp.term", "gpd.dens", "compfden.b1forintplaplace", "compfnum.vb1forintplaplace", "betahat") , envir = environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.claplace, seq.points, lys, betahat, tauvec,rhohat.l, xi.l,rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
     den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
    if(bdist == "Normal"){
      clusterExport(cl, c("seq.points", "lys", "XB.init", "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "compfnum.mb1forintp", "comp.term", "gpd.dens", "compfden.b1forintp", "compfnum.vb1forintp", "betahat"), envir =   environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.c, seq.points, lys, XB.init, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
     den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2, seq.points, lys,  betahat, tauvec,   rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
  }else{
    if(bdist == "Laplace"){
      cdf.out1 <-  t(sapply( 1:D,  ianum.mdp2.claplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.cdf <-  sapply(1:D,  iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
    if(bdist == "Normal"){
      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.cdf <-  sapply( 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    } 
  } 
  cdf.out <- cdf.out1/den.cdf
  bcond.ran <- replicate(Rb, gen.bcond(D,  cdf.out, seq.points))
  ######  sigma^2_b (1)
  mub.cond <- apply(bcond.ran, 1, mean)
  vebcond <-  apply(bcond.ran, 1, var)
  sig2bhat.update <- mean(bcond.ran^2)*D/(D-2)
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
   betahats.update <- estconstrained(dat.temp, as.formula("dev~X"),   tauvec)
  #betahats.update <- fit.update$coef
  XBhat.update1 <- cbind(1, lxN)%*%betahats.update
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}


updatebetasig2bfunSort <-  function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
  seq.points <- qnorm(tauvec, mean = 0, sd = sqrt(sig2bhat))
  lys <- dat.temp$Y
  rhohat.l <- gpdpar[1]; xi.l <- gpdpar[2]; rhohat.u <- gpdpar[3]; xi.u <- gpdpar[4]
  if(use.cl){
    if(bdist == "Laplace"){
      clusterExport(cl, c("seq.points", "lys", "XB.init", "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "compfnum.mb1forintplaplace", "comp.term", "gpd.dens", "compfden.b1forintplaplace", "compfnum.vb1forintplaplace", "betahat") , envir = environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.claplace, seq.points, lys, betahat, tauvec,rhohat.l, xi.l,rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
     den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
    if(bdist == "Normal"){
      clusterExport(cl, c("seq.points", "lys", "XB.init", "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "compfnum.mb1forintp", "comp.term", "gpd.dens", "compfden.b1forintp", "compfnum.vb1forintp", "betahat"), envir =   environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.c, seq.points, lys, XB.init, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
     den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2, seq.points, lys,  betahat, tauvec,   rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
  }else{
    if(bdist == "Laplace"){
      cdf.out1 <-  t(sapply( 1:D,  ianum.mdp2.claplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.cdf <-  sapply(1:D,  iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
    if(bdist == "Normal"){
      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.cdf <-  sapply( 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    } 
  } 
  cdf.out <- cdf.out1/den.cdf
  bcond.ran <- replicate(Rb, gen.bcond(D,  cdf.out, seq.points))
  ######  sigma^2_b (1)
  mub.cond <- apply(bcond.ran, 1, mean)
  vebcond <-  apply(bcond.ran, 1, var)
  sig2bhat.update <- mean(bcond.ran^2)*D/(D-2)
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
  fit.update <- rq(dev~X, data = dat.temp, tau = tauvec)
  betahats.update <- fit.update$coef
  XBhat.update <- cbind(1, lxN)%*%betahats.update
  XBhat.update1 <- t(apply( XBhat.update , 1,sort))
 
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}



updatebetasig2bfunConstrFix <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
  seq.points <- qnorm(tauvec, mean = 0, sd = sqrt(sig2bhat))
  lys <- dat.temp$Y
  rhohat.l <- gpdpar[1]; xi.l <- gpdpar[2]; rhohat.u <- gpdpar[3]; xi.u <- gpdpar[4]
  if(use.cl){
    if(bdist == "Laplace"){
      clusterExport(cl, c("seq.points", "lys", "XB.init", "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "compfnum.mb1forintplaplace", "comp.term", "gpd.dens", "compfden.b1forintplaplace", "compfnum.vb1forintplaplace", "betahat") , envir = environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.claplace, seq.points, lys, betahat, tauvec,rhohat.l, xi.l,rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
    if(bdist == "Normal"){
      clusterExport(cl, c("seq.points", "lys", "XB.init", "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "compfnum.mb1forintp", "comp.term", "gpd.dens", "compfden.b1forintp", "compfnum.vb1forintp", "betahat"), envir =   environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.c, seq.points, lys, XB.init, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2, seq.points, lys,  betahat, tauvec,   rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
  }else{
    if(bdist == "Laplace"){
      num.mub <-sapply( 1:D, FUN = ianum.mdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      num.vb <-sapply( 1:D, FUN = ianum.vdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      #      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
    if(bdist == "Normal"){
      num.mub <-sapply( 1:D, FUN = ianum.mdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      num.vb <-sapply( 1:D, FUN = ianum.vdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      #      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    } 
  } 
  mub.cond <- num.mub/den.mub
  vebcond <- num.vb/den.mub
  sig2bhat.update <-  mean(vebcond)*D/(D-(dim(cbind(1,lxN))[2]))
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
  betahats.update <- estconstrained(dat.temp, as.formula("dev~X"),   tauvec)
  #betahats.update <- fit.update$coef
  XBhat.update1 <- cbind(1, lxN)%*%betahats.update
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}




updatebetasig2bfunSortFix <-  function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
  seq.points <- qnorm(tauvec, mean = 0, sd = sqrt(sig2bhat))
  lys <- dat.temp$Y
  rhohat.l <- gpdpar[1]; xi.l <- gpdpar[2]; rhohat.u <- gpdpar[3]; xi.u <- gpdpar[4]
  if(use.cl){
    if(bdist == "Laplace"){
      clusterExport(cl, c("seq.points", "lys", "XB.init", "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "compfnum.mb1forintplaplace", "comp.term", "gpd.dens", "compfden.b1forintplaplace", "compfnum.vb1forintplaplace", "betahat") , envir = environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.claplace, seq.points, lys, betahat, tauvec,rhohat.l, xi.l,rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
    if(bdist == "Normal"){
      clusterExport(cl, c("seq.points", "lys", "XB.init", "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "compfnum.mb1forintp", "comp.term", "gpd.dens", "compfden.b1forintp", "compfnum.vb1forintp", "betahat"), envir =   environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.c, seq.points, lys, XB.init, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2, seq.points, lys,  betahat, tauvec,   rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
  }else{
    if(bdist == "Laplace"){
      num.mub <-sapply( 1:D, FUN = ianum.mdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      num.vb <-sapply( 1:D, FUN = ianum.vdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      #      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
    if(bdist == "Normal"){
      num.mub <-sapply( 1:D, FUN = ianum.mdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      num.vb <-sapply( 1:D, FUN = ianum.vdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      #      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    } 
  } 
  mub.cond <- num.mub/den.mub
  vebcond <- num.vb/den.mub
  sig2bhat.update <-  mean(vebcond)*D/(D-(dim(cbind(1,lxN))[2]))
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
  fit.update <- rq(dev~X, data = dat.temp, tau = tauvec)
  betahats.update <- fit.update$coef
  XBhat.update <- cbind(1, lxN)%*%betahats.update
  XBhat.update1 <- t(apply( XBhat.update , 1,sort))
  
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}




