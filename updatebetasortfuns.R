 

par.updatebetasig2bfunSortFix <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
  seq.points <- qnorm(tauvec, mean = 0, sd = sqrt(sig2bhat))
  lys <- dat.temp$Y
  rhohat.l <- gpdpar[1]; xi.l <- gpdpar[2]; rhohat.u <- gpdpar[3]; xi.u <- gpdpar[4]
    if(bdist == "Laplace"){
      clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XB.init", 
         "comp.term", "gpd.dens", "compfall.b1forintplaplace","lXs","dgpd") , envir = environment())
        out=parSapply(cl, 1:D, FUN = ianumden.mvdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
        }
    if(bdist == "Normal"){
      clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XB.init", 
        "comp.term", "gpd.dens", "compfall.b1forintp" ,"lXs","dgpd"), envir =   environment())
      out=parSapply(cl, 1:D, FUN = ianumden.mvdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)

    }

        num.mub=out[1,]
        num.vb=out[2,]
        den.mub=out[3,]
        
        mub.cond=num.mub/den.mub
        vebcond=num.vb/den.mub

  sig2bhat.update <-  mean(vebcond)*D/(D-(dim(cbind(1,lxN))[2]))
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
  fit.update <- rq(dev~X, data = dat.temp, tau = tauvec)
  betahats.update <- fit.update$coef
  XBhat.update <- cbind(1, lxN)%*%betahats.update
  XBhat.update1 <- t(apply( XBhat.update , 1,sort))
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}



 

par.updatebetasig2bfunIsoregFix <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
  seq.points <- qnorm(tauvec, mean = 0, sd = sqrt(sig2bhat))
  lys <- dat.temp$Y
  rhohat.l <- gpdpar[1]; xi.l <- gpdpar[2]; rhohat.u <- gpdpar[3]; xi.u <- gpdpar[4]
    if(bdist == "Laplace"){
      clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XB.init", 
         "comp.term", "gpd.dens", "compfall.b1forintplaplace","lXs","dgpd") , envir = environment())
        out=parSapply(cl, 1:D, FUN = ianumden.mvdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
        }
    if(bdist == "Normal"){
      clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XB.init", 
        "comp.term", "gpd.dens", "compfall.b1forintp" ,"lXs","dgpd"), envir =   environment())
      out=parSapply(cl, 1:D, FUN = ianumden.mvdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)

    }

        num.mub=out[1,]
        num.vb=out[2,]
        den.mub=out[3,]
        
        mub.cond=num.mub/den.mub
        vebcond=num.vb/den.mub

  sig2bhat.update <-  mean(vebcond)*D/(D-(dim(cbind(1,lxN))[2]))
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
  fit.update <- rq(dev~X, data = dat.temp, tau = tauvec)
  betahats.update <- fit.update$coef
  XBhat.update <- cbind(1, lxN)%*%betahats.update
  XBhat.update1 <- t(parApply(cl, X = XBhat.update , 1, function(x){ isoreg(tauvec, x)$yf }))
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}




