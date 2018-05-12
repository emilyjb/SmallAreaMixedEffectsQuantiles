updatebetasig2bfunfix <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
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
  sig2bhat.update <-  mean(vebcond)*D/(D-(dim(lxN)[2]+1))
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
  fit.update <- rq(dev~X, data = dat.temp, tau = tauvec)
  betahats.update <- fit.update$coef
  XBhat.update <- cbind(1, lxN)%*%betahats.update
  if(use.cl){
    XBhat.update1 <- t(parApply(cl, X = XBhat.update , 1, function(x){ isoreg(tauvec, x)$yf })) 
  }else{
    XBhat.update1 <- t(apply( XBhat.update , 1,  sort))
  }
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}


updatebetasig2bfunfixlog <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
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
      num.mub <-sapply( 1:D, FUN = ianum.mdp2log, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      num.vb <-sapply( 1:D, FUN = ianum.vdp2log, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
#      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2log, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    } 
  } 
  mub.cond <- num.mub/den.mub
  vebcond <- num.vb/den.mub
  sig2bhat.update <-  mean(vebcond)*D/(D-(dim(cbind(1,lxN))[2]))
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
  fit.update <- rq(dev~X, data = dat.temp, tau = tauvec)
  betahats.update <- fit.update$coef
  XBhat.update <- cbind(1, lxN)%*%betahats.update
  if(use.cl){
    XBhat.update1 <- t(parApply(cl, X = XBhat.update , 1, function(x){ isoreg(tauvec, x)$yf })) 
  }else{
    XBhat.update1 <- t(apply( XBhat.update , 1,  sort))
  }
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}

compfden.b1forintplog <- function(b1, Dpick, lyst, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        id <- which(areafac.pop[smc] == Dpick)
        pred.mat <- predfix[id,] + b1
        if(length(id) == 1){
          pred.mat <- matrix(pred.mat, nrow = 1)
        }
        if(length(id) == 1){
          get.tau.u <- min(which(lyst[id] <= pred.mat))
        }else{
          get.tau.u <- apply(lyst[id] <= pred.mat, 1, function(x){ min(which(x))})
        }
        ll.vec <- sapply(1:length(id), comp.term, pred.mat,get.tau.u, lyst[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
       ll <- sum(log(ll.vec)) + log(dnorm(b1, mean = 0, sd = sdb))
	if(exp(ll)==0){
		exp(ll - log(0.5^length(id)))
	}else{
		exp(ll)
	}
}

compfnum.b1forintplog <- function(b1, Dpick, lyst, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        id <- which(areafac.pop[smc] == Dpick)
        pred.mat <- predfix[id,] + b1
        if(length(id) == 1){
          pred.mat <- matrix(pred.mat, nrow = 1)
        }
        if(length(id) == 1){
          get.tau.u <- min(which(lyst[id] <= pred.mat))
        }else{
          get.tau.u <- apply(lyst[id] <= pred.mat, 1, function(x){ min(which(x))})
        }
        ll.vec <- sapply(1:length(id), comp.term, pred.mat,get.tau.u, lyst[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
       ll1 <- sum(log(ll.vec)) + log(dnorm(b1, mean = 0, sd = sdb)) 
	if(exp(ll1) == 0){ 
		 exp(ll1  - log(0.5^length(id)))*b1
	}else{
		exp(ll1 )*b1
	}
}


compfnum.vb1forintplog <- function(b1, Dpick, lyst, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        id <- which(areafac.pop[smc] == Dpick)
        pred.mat <- predfix[id,] + b1
        if(length(id) == 1){
          pred.mat <- matrix(pred.mat, nrow = 1)
        }
        if(length(id) == 1){
          get.tau.u <- min(which(lyst[id] <= pred.mat))
        }else{
          get.tau.u <- apply(lyst[id] <= pred.mat, 1, function(x){ min(which(x))})
        }
        ll.vec <- sapply(1:length(id), comp.term, pred.mat,get.tau.u, lyst[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
       ll1 <- sum(log(ll.vec)) + log(dnorm(b1, mean = 0, sd = sdb))  
	 if(exp(ll1) == 0){
		 exp(ll1  - log(0.5^length(id)))*b1^2
	}else{
		exp(ll1)*b1^2
	}	
}


# compfnum.vb1forintp <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
#        id <- which(areafac.pop[smc] == Dpick)
#        pred.mat <- predfix[id,] + b1
#         if(length(id) == 1){
#          pred.mat <- matrix(pred.mat, nrow = 1)
#        }
#        if(length(id) == 1){
#          get.tau.u <- min(which(lyst[id] <= pred.mat))
#        }else{
#          get.tau.u <- apply(lyst[id] <= pred.mat, 1, function(x){ min(which(x))})
#        }
#        ll.vec <- sapply(1:length(lys[id]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
#        prod(ll.vec)*dnorm(b1, mean = 0, sd = sdb)*b1^2
#}


iaden.mdp2log <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        vec1 <- sapply(seq.points, compfden.b1forintplog, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix) 
  vec2 <- c(0, vec1[-length(vec1)])
  midpts <- (vec1 + vec2)/2
  widths <- diff(seq.points)
  sum(widths*(midpts[-1]))

}

ianum.mdp2log <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        vec1 <- sapply(seq.points, compfnum.b1forintplog, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix) 
  vec2 <- c(0, vec1[-length(vec1)])
  midpts <- (vec1 + vec2)/2
  widths <- diff(seq.points)
  sum(widths*(midpts[-1]))

}


ianum.vdp2log <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        vec1 <- sapply(seq.points, compfnum.vb1forintplog, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix) 
  vec2 <- c(0, vec1[-length(vec1)])
  midpts <- (vec1 + vec2)/2
  widths <- diff(seq.points)
  sum(widths*(midpts[-1]))

}


 

updatebetasig2bfunfixcheckV0 <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
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
##### cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    }
    if(bdist == "Normal"){
      num.mub <-sapply( 1:D, FUN = ianum.mdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      num.vb <-sapply( 1:D, FUN = ianum.vdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
##### cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    } 
  }
   
  mub.cond <- num.mub/den.mub
  vebcond <- num.vb/den.mub
  sig2bhat.update <-  mean(vebcond)*D/(D-(length(varspick)+1))
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
  fit.update <- rq(dev~X, data = dat.temp, tau = tauvec)
  betahats.update <- fit.update$coef
  XBhat.update <- cbind(1, lxN)%*%betahats.update
  if(use.cl){
    XBhat.update1 <- t(parApply(cl, X = XBhat.update , 1, function(x){ isoreg(tauvec, x)$yf })) 
  }else{
    XBhat.update1 <- t(apply( XBhat.update , 1,  sort))
  }
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}




 predJWTransFix <-  function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XBhat.update1, Gs, lxN, use.cl= TRUE, Rb = 150 , trunc = FALSE, lam){
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
     num.mub <-sapply( 1:D, FUN = ianum.mdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
      num.vb <-sapply( 1:D, FUN = ianum.vdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
    if(bdist == "Normal"){
     num.mub <-sapply( 1:D, FUN = ianum.mdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
      num.vb <-sapply( 1:D, FUN = ianum.vdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    } 
  } 
  mub.cond <- num.mub/den.mub
  vebcond <- num.vb/den.mub
    
  all.q.J <- XBhat.update1 + as.vector(GN%*%mub.cond)
  t.min.J <- min(all.q.J)
  t.max.J <- max(all.q.J)
  t.seq.J <- seq(t.min.J, t.max.J, length = 100)
  
  if(trunc){
    
    all.q.J[all.q.J <= 0] <- 0
    
  }
  
  all.q.Jtrans <- t(apply(all.q.J, 1, hinvtransG1, lam) - 0.001)
 # all.q.Jtrans[smc,] <- lys
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


 predJWFix <-  function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XBhat.update1, Gs, lxN, use.cl= TRUE, Rb = 150 , trunc = FALSE, lam){
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
     num.mub <-sapply( 1:D, FUN = ianum.mdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
      num.vb <-sapply( 1:D, FUN = ianum.vdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
    if(bdist == "Normal"){
     num.mub <-sapply( 1:D, FUN = ianum.mdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
      num.vb <-sapply( 1:D, FUN = ianum.vdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    } 
  } 
  mub.cond <- num.mub/den.mub
  vebcond <- num.vb/den.mub
    
  all.q.J <- XBhat.update1 + as.vector(GN%*%mub.cond)
  t.min.J <- min(all.q.J)
  t.max.J <- max(all.q.J)
  t.seq.J <- seq(t.min.J, t.max.J, length = 100)
  
  if(trunc){
    
    all.q.J[all.q.J <= 0] <- 0
    
  }
  
  all.q.Jtrans <-  all.q.J 
#  all.q.Jtrans[smc,] <- lys
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



 predJWFixlog <-  function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XBhat.update1, Gs, lxN, use.cl= TRUE, Rb = 150 , trunc = FALSE, lam){
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
     num.mub <-sapply( 1:D, FUN = ianum.mdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
      num.vb <-sapply( 1:D, FUN = ianum.vdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
    if(bdist == "Normal"){
     num.mub <-sapply( 1:D, FUN = ianum.mdp2log, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
      num.vb <-sapply( 1:D, FUN = ianum.vdp2log, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2log, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    } 
  } 
  mub.cond <- num.mub/den.mub
  vebcond <- num.vb/den.mub
    
  all.q.J <- XBhat.update1 + as.vector(GN%*%mub.cond)
  t.min.J <- min(all.q.J)
  t.max.J <- max(all.q.J)
  t.seq.J <- seq(t.min.J, t.max.J, length = 100)
  
  if(trunc){
    
    all.q.J[all.q.J <= 0] <- 0
    
  }
  
  all.q.Jtrans <-  all.q.J 
 # all.q.Jtrans[smc,] <- lys
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


updatebetasig2bfunfixsort <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
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
  sig2bhat.update <-  mean(vebcond)*D/(D-dim(cbind(1,lxN))[2])
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
  fit.update <- rq(dev~X, data = dat.temp, tau = tauvec)
  betahats.update <- fit.update$coef
  XBhat.update <- cbind(1, lxN)%*%betahats.update
   XBhat.update1 <- t(apply(XBhat.update,1, sort))
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}




updatebetasig2bfunfixsortmatch <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
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
      num.mub <-sapply( 1:D, FUN = ianum.mdpfixMatch, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      num.vb <-sapply( 1:D, FUN = ianum.vdpfixMatch, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init)
      #      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XB.init))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2Match, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
    } 
  } 
  mub.cond <- num.mub/den.mub
  vebcond <- num.vb/den.mub
  sig2bhat.update <-  mean(vebcond)*D/(D-dim(cbind(1,lxN))[2])
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
  fit.update <- rq(dev~X, data = dat.temp, tau = tauvec)
  betahats.update <- fit.update$coef
  XBhat.update <- cbind(1, lxN)%*%betahats.update
  XBhat.update1 <- t(apply(XBhat.update,1, sort))
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}


predJWTransFixB0 <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XBhat.update1, Gs, lxN, use.cl= TRUE, Rb = 150 , trunc = FALSE, lam, B = 0){
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
     num.mub <-sapply( 1:D, FUN = ianum.mdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
      num.vb <-sapply( 1:D, FUN = ianum.vdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
    if(bdist == "Normal"){
     num.mub <-sapply( 1:D, FUN = ianum.mdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
      num.vb <-sapply( 1:D, FUN = ianum.vdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    } 
  } 
  mub.cond <- num.mub/den.mub
  vebcond <- num.vb/den.mub
    
  all.q.J <- XBhat.update1 + as.vector(GN%*%mub.cond)
  t.min.J <- min(all.q.J)
  t.max.J <- max(all.q.J)
  t.seq.J <- seq(t.min.J, t.max.J, length = 100)
  
  if(trunc){
    
    all.q.J[all.q.J <= 0] <- 0
    
  }
  
  all.q.Jtrans <- t(apply(all.q.J, 1, hinvtransG1, lam) - B)
#  all.q.Jtrans[smc,] <- lys
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





predJWTransFixB0ExpOld <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XBhat.update1, Gs, lxN, use.cl= TRUE, Rb = 150 , trunc = FALSE, lam, B = 0){
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
     num.mub <-sapply( 1:D, FUN = ianum.mdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
      num.vb <-sapply( 1:D, FUN = ianum.vdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
    if(bdist == "Normal"){
     num.mub <-sapply( 1:D, FUN = ianum.mdpfixexp, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
      num.vb <-sapply( 1:D, FUN = ianum.vdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    } 
  } 
  mub.cond <- num.mub/den.mub
  vebcond <- num.vb/den.mub
    
  all.q.J <- all.q.Jtrans <- exp(XBhat.update1)*as.vector(GN%*%mub.cond)
  t.min.J <- min(all.q.J)
  t.max.J <- max(all.q.J)
  t.seq.J <- seq(t.min.J, t.max.J, length = 100)
  
  if(trunc){
    
    all.q.J[all.q.J <= 0] <- 0
    
  }
  
#  all.q.Jtrans <- t(apply(all.q.J, 1, hinvtransG1, lam) - B)
#  all.q.Jtrans[smc,] <- lys
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




