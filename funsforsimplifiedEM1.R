intparfunfordatmulticov <- function(dat.temp, X, areafac.samp,D, Gs, tauvec, lxN,  use.cl = TRUE){
  
  mm <- model.matrix(~lxN[smc,] + as.factor(areafac.pop[smc]))
  mmsolve <- NULL
  try(mmsolve <- solve(t(mm)%*%mm))
  if(!is.null(mmsolve)){        
     lm1med <- rq(dat.temp$Y~lxN[smc,] + as.factor(areafac.pop[smc]), tau = 0.5)
     bhats <-  summary.rq(lm1med, se = "ker")$coefficients[,"Value"][-c(1:(dim(X)[2]+1))]
     vbhat <-  summary.rq(lm1med, se = "ker")$coefficients[,"Std. Error"][-c(1:(dim(X)[2]+1))]^2
  }else{
      lm1med <- rq(dat.temp$Y~ as.factor(areafac.pop[smc]), tau = 0.5)
      bhats <-  summary.rq(lm1med, se = "ker")$coefficients[,"Value"][-c(1)]
      vbhat <-  summary.rq(lm1med, se = "ker")$coefficients[,"Std. Error"][-c(1)]^2
  }       

  sig2bhat <- betahat.glsfun(vbhat, bhats , 2, matrix(1,nrow = D-1), bhats)[[2]]

  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%c( bhats,  -sum(bhats)))
  fit.init <- rq(dat.temp$dev~X,  tau = tauvec)
  paravec <- c(as.vector(fit.init$coef), c( bhats,  -sum(bhats)))
  XBhat.init <- cbind(1, lxN)%*%fit.init$coef
  if(use.cl){
    clusterExport(cl, "tauvec")
    XBhat.init1 <-  t( parApply(cl, XBhat.init, 1, function(x){ isoreg(tauvec, x)$yf }))
  }else{
    XBhat.init1 <- t( apply( XBhat.init, 1, function(x){ isoreg(tauvec, x)$yf }))
  }
  list(beta = fit.init$coef, sig2bhat= sig2bhat, bhats = bhats, XBinit = XBhat.init1)
}



intparfunfordat <- function(dat.temp, X, areafac.samp,D, Gs, tauvec, lxN, use.cl = TRUE){
  
  mm <- model.matrix(~lxN[smc] + as.factor(areafac.pop[smc]))
  mmsolve <- NULL
  try(mmsolve <- solve(t(mm)%*%mm))
  if(!is.null(mmsolve)){	
     lm1med <- rq(dat.temp$Y~lxN[smc] + as.factor(areafac.pop[smc]), tau = 0.5)
     bhats <-  summary.rq(lm1med, se = "ker")$coefficients[,"Value"][-c(1,2)]
     vbhat <-  summary.rq(lm1med, se = "ker")$coefficients[,"Std. Error"][-c(1,2)]^2
  }else{
      lm1med <- rq(dat.temp$Y~ as.factor(areafac.pop[smc]), tau = 0.5)
		bhats <-  summary.rq(lm1med, se = "ker")$coefficients[,"Value"][-c(1)]
		vbhat <-  summary.rq(lm1med, se = "ker")$coefficients[,"Std. Error"][-c(1)]^2
	}	

  sig2bhat <- betahat.glsfun(vbhat, bhats , 2, matrix(1,nrow = D-1), bhats)[[2]]

  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%c( bhats,  -sum(bhats)))
  fit.init <- rq(dev~X, data = dat.temp, tau = tauvec)
  paravec <- c(as.vector(fit.init$coef), c( bhats,  -sum(bhats)))
  XBhat.init <- cbind(1, lxN)%*%fit.init$coef
  if(use.cl){
    clusterExport(cl, "tauvec")
    XBhat.init1 <-  t( parApply(cl, XBhat.init, 1, function(x){ isoreg(tauvec, x)$yf }))
  }else{
    XBhat.init1 <- t( apply( XBhat.init, 1, function(x){ isoreg(tauvec, x)$yf }))
  }
  list(beta = fit.init$coef, sig2bhat= sig2bhat, bhats = bhats, XBinit = XBhat.init1)
}



intparfun <- function(dat.temp, X, areafac.samp,D, Gs, tauvec, lxN, use.cl = TRUE){
  
  lm1med <- rq(dat.temp$Y~X + as.factor(areafac.samp), tau = 0.5)
  bhats <-  summary.rq(lm1med, se = "ker")$coefficients[,"Value"][-c(1,2)]
  vbhat <-  summary.rq(lm1med, se = "ker")$coefficients[,"Std. Error"][-c(1,2)]^2
  sig2bhat <- betahat.glsfun(vbhat, bhats , 2, matrix(1,nrow = D-1), bhats)[[2]]

  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%c( bhats,  -sum(bhats)))
  fit.init <- rq(dev~X, data = dat.temp, tau = tauvec)
  paravec <- c(as.vector(fit.init$coef), c( bhats,  -sum(bhats)))
  XBhat.init <- cbind(1, lxN)%*%fit.init$coef
  if(use.cl){
    clusterExport(cl, "tauvec")
    XBhat.init1 <-  t( parApply(cl, XBhat.init, 1, function(x){ isoreg(tauvec, x)$yf }))
  }else{
    XBhat.init1 <- t( apply( XBhat.init, 1, function(x){ isoreg(tauvec, x)$yf }))
  }
  list(beta = fit.init$coef, sig2bhat= sig2bhat, bhats = bhats, XBinit = XBhat.init1)
}

#### Note: XB.init is XB matrix for sampled units
gpdparMOM <- function(XB.init, dat.temp, tauvec){
  ### Estimate GPD parameters
  dev.lower <- XB.init[,2] - XB.init[,1]
  dev.upper <- XB.init[,length(tauvec)] - XB.init[,length(tauvec)-1]
  
  #rhohat.l <- max(c((tauvec[1] + tauvec[2])/(2*nsamp*(tauvec[2] - tauvec[1]))*sum(dev.lower),10^-10))
  #rhohat.u <- max(c((2-(tauvec[length(tauvec)] + tauvec[length(tauvec)-1]))/(2*nsamp*(tauvec[length(tauvec)] - tauvec[length(tauvec)-1]))*sum(dev.upper), 10^-10))
  
  lij <- (XB.init[,2] + XB.init[,1])/2
  uij <- (XB.init[,length(tauvec)] + XB.init[,length(tauvec)-1])/2
  
  ex.l <- dat.temp$Y[dat.temp$Y <= lij]
  ex.u <- dat.temp$Y[dat.temp$Y >= uij]
  
  du <- ex.u - uij[dat.temp$Y >= uij]
  if(length(du) <= 2){
    du <- dat.temp$Y[dat.temp$Y >= sort(dat.temp$Y)[length(dat.temp$Y - 1)]]
  }
  xi.u <- (1 - mean(du)^2/var(du))/2
  rhohat.u <- (1 - xi.u)*mean(du)
  
  dl <- lij[dat.temp$Y <= lij] - ex.l
  if(length(dl) <= 1){
    dl <- dat.temp$Y[dat.temp$Y <= sort(dat.temp$Y)[2]]
  }
  xi.l <- (1 - mean(dl)^2/var(dl))/2
  rhohat.l <- (1 - xi.l)*mean(dl)
  
  c( rhohat.l,xi.l,rhohat.u,  xi.u)

}

gpdparMatch <- function(XB.init, dat.temp, tauvec){
	nsamp <- nrow(XB.init)
	dev.lower <- XB.init[,2] - XB.init[,1]
	dev.upper <- XB.init[,length(tauvec)] - XB.init[,length(tauvec)-1]

	rhohat.l <- max(c((tauvec[1] + tauvec[2])/(2*nsamp*(tauvec[2] - tauvec[1]))*sum(dev.lower),10^-10))
	rhohat.u <- max(c((2-(tauvec[length(tauvec)] + tauvec[length(tauvec)-1]))/(2*nsamp*(tauvec[length(tauvec)] - tauvec[length(tauvec)-1]))*sum(dev.upper), 10^-10))

	lij <- (XB.init[,2] + XB.init[,1])/2
	uij <- (XB.init[,length(tauvec)] + XB.init[,length(tauvec)-1])/2

	ex.l <- dat.temp$Y[dat.temp$Y <= lij]
	ex.u <- dat.temp$Y[dat.temp$Y >= uij]

	du <- ex.u - uij[dat.temp$Y >= uij]
	if(length(du) <= 2){
		du <- dat.temp$Y[dat.temp$Y >= sort(dat.temp$Y)[length(dat.temp$Y - 1)]]
	}

	xifit.u <- optimize(gpdxi, interval = c(0, 100), maximum = TRUE, dvec = du, rho = rhohat.u, taul = tauvec[length(tauvec)-1], tauu = tauvec[length(tauvec)])
	xi.u <- xifit.u$maximum

	dl <- lij[dat.temp$Y <= lij] - ex.l
	if(length(dl) <= 1){
		dl <- dat.temp$Y[dat.temp$Y <= sort(dat.temp$Y)[2]]
	}

	xifit.l <- optimize(gpdxi, interval = c(0, 100), maximum = TRUE, dvec = dl, rho = rhohat.l, taul = tauvec[1], tauu = tauvec[2])
	xi.l <- xifit.l$maximum

      c( rhohat.l,xi.l,  rhohat.u, xi.u)

}



### Note: XB.init is the XB matrix for only the sampled observations
updatebetasig2bfun <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
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
  if(use.cl){
    XBhat.update1 <- t(parApply(cl, X = XBhat.update , 1, function(x){ isoreg(tauvec, x)$yf })) 
  }else{
    XBhat.update1 <- t(apply( XBhat.update , 1, function(x){ isoreg(tauvec, x)$yf }))
  }
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}


predJW <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XBhat.update1, Gs, lxN, use.cl= TRUE, Rb = 150 , trunc = FALSE){
  seq.points <- qnorm(tauvec, mean = 0, sd = sqrt(sig2bhat))
  lys <- dat.temp$Y
  rhohat.l <- gpdpar[1]; xi.l <- gpdpar[2]; rhohat.u <- gpdpar[3]; xi.u <- gpdpar[4]
  if(use.cl){
    if(bdist == "Laplace"){
      clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XBhat.update1", "compfnum.mb1forintplaplace", "comp.term", "gpd.dens", "compfden.b1forintplaplace", "compfnum.vb1forintplaplace") , envir = environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.claplace, seq.points, lys, betahat, tauvec,rhohat.l, xi.l,rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
    if(bdist == "Normal"){
      clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XBhat.update1", "compfnum.mb1forintp", "comp.term", "gpd.dens", "compfden.b1forintp", "compfnum.vb1forintp"), envir =   environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.c, seq.points, lys,betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
  }else{
    if(bdist == "Laplace"){
      cdf.out1 <- t(sapply( 1:D,  ianum.mdp2.claplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.cdf <- sapply(1:D,  iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
    if(bdist == "Normal"){
      cdf.out1 <- t(sapply( 1:D,  ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.cdf <- sapply( 1:D,   iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    } 
  } 
  cdf.out <- cdf.out1/den.cdf
  bcond.ran <- replicate(Rb, gen.bcond(D,  cdf.out, seq.points))
  mub.cond <- apply(bcond.ran, 1, mean)
  vebcond <- apply(bcond.ran, 1, var)
  
  all.q.J <- XBhat.update1 + as.vector(GN%*%mub.cond)
  t.min.J <- min(all.q.J)
  t.max.J <- max(all.q.J)
  t.seq.J <- seq(t.min.J, t.max.J, length = 100)
  
  if(trunc){
	
 	all.q.J[all.q.J <= 0] <- 0

 }

	
  w25J <- sapply(1:D, function(i){ quantile(as.vector( all.q.J[areafac.pop == i,]), prob = 0.25)})
  w50J <-  sapply(1:D, function(i){ quantile(as.vector( all.q.J[areafac.pop == i,]), prob = 0.5)})
  w75J <-  sapply(1:D, function(i){ quantile(as.vector( all.q.J[areafac.pop == i,]), prob = 0.75)})
  w90J <- sapply(1:D, function(i){ quantile(as.vector( all.q.J[areafac.pop == i,]), prob = 0.90)})
  w10J <-  sapply(1:D, function(i){ quantile(as.vector( all.q.J[areafac.pop == i,]), prob = 0.10)})
  
  meanJ <- sapply(1:D, function(i){ mean(as.vector( all.q.J[areafac.pop == i,]))})
  varJ <-  sapply(1:D, function(i){ var(as.vector( all.q.J[areafac.pop == i,]))})
  
  mubJ <-  mub.cond
  vebJ <-  vebcond
  
   list( w25J = w25J, w50J = w50J, w75J = w75J, w90J = w90J, w10J = w10J, meanJ = meanJ, varJ = varJ,mubJ = mubJ, vebJ = vebJ, all.q.J = all.q.J)
   
}



predJWTrans <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XBhat.update1, Gs, lxN, use.cl= TRUE, Rb = 150 , trunc = FALSE, lam){
  seq.points <- qnorm(tauvec, mean = 0, sd = sqrt(sig2bhat))
  lys <- dat.temp$Y
  rhohat.l <- gpdpar[1]; xi.l <- gpdpar[2]; rhohat.u <- gpdpar[3]; xi.u <- gpdpar[4]
  if(use.cl){
    if(bdist == "Laplace"){
      clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XBhat.update1", "compfnum.mb1forintplaplace", "comp.term", "gpd.dens", "compfden.b1forintplaplace", "compfnum.vb1forintplaplace") , envir = environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.claplace, seq.points, lys, betahat, tauvec,rhohat.l, xi.l,rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
    if(bdist == "Normal"){
      clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XBhat.update1", "compfnum.mb1forintp", "comp.term", "gpd.dens", "compfden.b1forintp", "compfnum.vb1forintp"), envir =   environment())
      cdf.out1 <- t(parSapply(cl, X = 1:D, FUN = ianum.mdp2.c, seq.points, lys,betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.cdf <- parSapply(cl, X = 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
  }else{
    if(bdist == "Laplace"){
      cdf.out1 <- t(sapply( 1:D,  ianum.mdp2.claplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.cdf <- sapply(1:D,  iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
    if(bdist == "Normal"){
      cdf.out1 <- t(sapply( 1:D,  ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.cdf <- sapply( 1:D,   iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    } 
  } 
  cdf.out <- cdf.out1/den.cdf
  bcond.ran <- replicate(Rb, gen.bcond(D,  cdf.out, seq.points))
  mub.cond <- apply(bcond.ran, 1, mean)
  vebcond <- apply(bcond.ran, 1, var)
  
  all.q.J <- XBhat.update1 + as.vector(GN%*%mub.cond)
  t.min.J <- min(all.q.J)
  t.max.J <- max(all.q.J)
  t.seq.J <- seq(t.min.J, t.max.J, length = 100)
  
  if(trunc){
    
    all.q.J[all.q.J <= 0] <- 0
    
  }
  
  all.q.Jtrans <- t(apply(all.q.J, 1, hinvtransG1, lam) - 0.001)
  w25J <- sapply(1:D, function(i){ quantile(as.vector( all.q.Jtrans[areafac.pop == i,]), prob = 0.25)})
  w50J <-  sapply(1:D, function(i){ quantile(as.vector( all.q.Jtrans[areafac.pop == i,]), prob = 0.5)})
  w75J <-  sapply(1:D, function(i){ quantile(as.vector( all.q.Jtrans[areafac.pop == i,]), prob = 0.75)})
  w90J <- sapply(1:D, function(i){ quantile(as.vector( all.q.Jtrans[areafac.pop == i,]), prob = 0.90)})
  w10J <-  sapply(1:D, function(i){ quantile(as.vector( all.q.Jtrans[areafac.pop == i,]), prob = 0.10)})
  
  meanJ <- sapply(1:D, function(i){ mean(as.vector( all.q.Jtrans[areafac.pop == i,]))})
  varJ <-  sapply(1:D, function(i){ var(as.vector( all.q.Jtrans[areafac.pop == i,]))})
  
  mubJ <-  mub.cond
  vebJ <-  vebcond
  
  list( w25J = w25J, w50J = w50J, w75J = w75J, w90J = w90J, w10J = w10J, meanJ = meanJ, varJ = varJ,mubJ = mubJ, vebJ = vebJ, all.q.J = all.q.J, all.q.Jtrans = all.q.Jtrans)
  
}




  


