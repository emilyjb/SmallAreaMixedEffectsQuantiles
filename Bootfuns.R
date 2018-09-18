bootgensimple3 <- function(sig2bhat.update, GN, XBhat.update1, CNis, nis, D, areafac.pop, rholxilrhouxiu, smc, b.dist){
if(b.dist == "Normal"){
  ahats.b <- rnorm(D, mean = 0, sd = sqrt(sig2bhat.update ))
}else{
  Ul <- runif(D, -0.5, 0.5)
  b <- sqrt(sig2bhat.update/2)
  ahats.b <- -b*sign(Ul)*log(1 - 2*abs(Ul))
}
  all.y.b <- XBhat.update1 + as.vector(GN%*%ahats.b)
  lyN.b <- areapredraninterp(all.y.b, rholxilrhouxiu, tauvec)	
  q25boot  <-  tapply(lyN.b, areafac.pop, quantile, prob = 0.25)
  q50boot <-   tapply(lyN.b, areafac.pop, quantile, prob = 0.5)
  q75boot <-  tapply(lyN.b, areafac.pop, quantile, prob = 0.75)

  lys.b <- lyN.b[smc]

  list(lys.b, smc, q25boot = q25boot, q50boot = q50boot, q75boot = q75boot)

}


areapredraninterp  <- function(all.q.Jpop, rholxilrhouxiu, tauvec){
	tausRan <- runif(nrow(all.q.Jpop))
	lowerLimTau <- mean(tauvec[c(1,2)])
	upperLimTau <- mean(tauvec[c(length(tauvec),length(tauvec)-1)])
	tausRanLower <- tausRan[tausRan < lowerLimTau]/lowerLimTau
	tausRanUpper <- (tausRan[tausRan > upperLimTau]  - (tauvec[length(tauvec)]+tauvec[length(tauvec)-1])/2)/(1 - (tauvec[length(tauvec)]+tauvec[length(tauvec)-1])/2)
	ystarlower <- -qgpd(tausRanLower, rholxilrhouxiu[1], rholxilrhouxiu[2]) + apply(all.q.Jpop[which(tausRan < lowerLimTau),c(1,2)], 1, mean)
	ystarupper <-  qgpd(tausRanUpper, rholxilrhouxiu[3], rholxilrhouxiu[4]) + apply(all.q.Jpop[which(tausRan >upperLimTau),c(length(tauvec)-1,length(tauvec))], 1, mean)
	ks <- sapply(tausRan, function(t){ 
	out=c()
	if (sum(t<=tauvec)==0) out=Inf
      if (sum(t<=tauvec)>0) out=min(which(t <= tauvec))
	return(out)})
	ystar1 <- sapply(1:nrow(all.q.Jpop), interpkreverse,  ks, tausRan, tauvec,all.q.Jpop)
	ystar1[ tausRan < lowerLimTau] <- ystarlower
	ystar1[ tausRan > upperLimTau] <- ystarupper
	ystar1
}


interpk <- function(  i,ks, ys,  all.q.Jsamp){
	if(ks[i] == Inf){
		return(1)
	}
	if(ks[i] == 1){
		return(0)
	}
      checkdiff <- all.q.Jsamp[i,ks[i]] - all.q.Jsamp[i,ks[i]-1]
	if(abs(checkdiff) < 10^-5){
		all.q.Jsamp[i,ks[i]-1]
	}else{
		tauvec[ks[i]-1] + (ys[i] - all.q.Jsamp[i,ks[i]-1])*( tauvec[ks[i]] - tauvec[ks[i]-1])/(all.q.Jsamp[i,ks[i]] - all.q.Jsamp[i,ks[i]-1])
	}
}


interpkreverse <- function(i, ks, taus, tauvec, all.q.Jpop){
	if(ks[i] == Inf){
		return(1)
	}
	if(ks[i] == 1){
		return(0)
	}
	checkdiff <- all.q.Jpop[i,ks[i]] - all.q.Jpop[i,ks[i]-1]
	if(abs(checkdiff) < 10^-5){
		all.q.Jpop[i,ks[i]-1]
	}else{
		all.q.Jpop[i,ks[i]-1] + (taus[i] - tauvec[ks[i]-1])*( (all.q.Jpop[i,ks[i]] - all.q.Jpop[i,ks[i]-1])/( tauvec[ks[i]] - tauvec[ks[i]-1]))
	}
}


par.bootestsimple1=function(dat.tempb, X, areafac.samp, D, Gs, tauvec,  Rb, lxN, b.dist, areafac.pop, smc, use.cl, trunc =FALSE){

  initpars <- intparfun(dat.tempb, X, areafac.samp,D, Gs,tauvec,  lxN, use.cl)

  rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.tempb, tauvec)

  XBbetasiguupdate <- par.updatebetasig2bfun(tauvec,dat.tempb, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl, Rb )

  sig2bhatupdate <- XBbetasiguupdate[[3]]
  betahat <- XBbetasiguupdate[[2]]
  XBhatupdate <- XBbetasiguupdate[[1]]

  XBbetasiguupdate2 <- par.updatebetasig2bfun(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat , rholxilrhouxiu , areafac.pop, smc, XBhatupdate[smc,], Gs,lxN, use.cl, Rb = 150 )

  sig2bhatupdate <- XBbetasiguupdate2[[3]]
  betahat <- XBbetasiguupdate2[[2]]
  XBhatupdate <- XBbetasiguupdate2[[1]]

  areapred <- par.predJW(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl, Rb = 150, trunc )

  qhat25b <- areapred[[1]]
  qhat50b <- areapred[[2]]
  qhat75b <- areapred[[3]]
  
  list(qhat25b, qhat50b, qhat75b, sig2bhatupdate, betahat, XBhatupdate)
  
}


par.bootestsimple1Rev =function(dat.tempb, X, areafac.samp, D, Gs, tauvec,  Rb, lxN, b.dist, areafac.pop, smc, use.cl, trunc =FALSE){

  initpars <- intparfun(dat.tempb, X, areafac.samp,D, Gs,tauvec,  lxN, use.cl)

  rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.tempb, tauvec)

  XBbetasiguupdate <- par.updatebetasig2bfunSortFix(tauvec,dat.tempb, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl, Rb )

  sig2bhatupdate <- XBbetasiguupdate[[3]]
  betahat <- XBbetasiguupdate[[2]]
  XBhatupdate <- XBbetasiguupdate[[1]]

  rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.tempb, tauvec)

  XBbetasiguupdate <- par.updatebetasig2bfunSortFix(tauvec,dat.tempb, sig2bhatupdate, b.dist, betahat,rholxilrhouxiu , areafac.pop, smc, XBhatupdate[smc,], Gs, lxN, use.cl, Rb )

  sig2bhatupdate <- XBbetasiguupdate[[3]]
  betahat <- XBbetasiguupdate[[2]]
  XBhatupdate <- XBbetasiguupdate[[1]]

  rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.tempb, tauvec)


  areapred <- par.predJWFix(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl, Rb = 150, trunc )

  qhat25b <- areapred[[1]]
  qhat50b <- areapred[[2]]
  qhat75b <- areapred[[3]]
  
  list(qhat25b, qhat50b, qhat75b, sig2bhatupdate, betahat, XBhatupdate)
  
}


par.updatebetasig2bfun=function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
  seq.points <- qnorm(tauvec, mean = 0, sd = sqrt(sig2bhat))
  lys <- dat.temp$Y
  rhohat.l <- gpdpar[1]; xi.l <- gpdpar[2]; rhohat.u <- gpdpar[3]; xi.u <- gpdpar[4]
    if(bdist == "Laplace"){
	  clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XB.init", 
	 "comp.term", "gpd.dens", "compfden.b1forintplaplace","lXs","dgpd") , envir = environment())
    	cdf.out=t(parSapply(cl, 1:D, FUN = ianumden.mdp2.claplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init))
       }
    if(bdist == "Normal"){
   clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat", "XB.init", 
	"comp.term", "gpd.dens", "compfden.b1forintp" ,"lXs","dgpd"), envir =   environment())
      cdf.out=t(parSapply(cl, 1:D, FUN = ianumden.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init))
    } 
  bcond.ran <- replicate(Rb, gen.bcond(D,  cdf.out, seq.points))
  ######  sigma^2_b (1)
  mub.cond <- apply(bcond.ran, 1, mean)
  vebcond <-  apply(bcond.ran, 1, var)
  sig2bhat.update <- mean(bcond.ran^2)*D/(D-2)
  dat.temp$dev <- dat.temp$Y - as.vector(Gs%*%mub.cond)
  fit.update <- rq(dev~X, data = dat.temp, tau = tauvec)
  betahats.update <- fit.update$coef
  XBhat.update <- cbind(1, lxN)%*%betahats.update
  XBhat.update1 <- t(apply( XBhat.update , 1, function(x){ isoreg(tauvec, x)$yf }))
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}


compfden.b1forintplaplace <- function(b1, Dpick, lyst, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- predfix[id,] + b1
	get.tau.u <- apply(lyst[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(id), comp.term, pred.mat,get.tau.u, lyst[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dlp(b1, mean = 0, sd = sdb)
}


ianumden.mdp2.c=function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        vec1 <- sapply(seq.points, compfden.b1forintp, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)
       vec2 <- c(0, vec1[-length(vec1)])
    midpts <- (vec1 + vec2)/2
   widths <- diff(seq.points)
  out1=cumsum(widths*(midpts[-1]))
  out2=sum(widths*(midpts[-1]))
	return(out1/out2)
}


ianumden.mdp2.claplace=function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        vec1 <- sapply(seq.points, compfden.b1forintplaplace, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)
          vec2 <- c(0, vec1[-length(vec1)])
          midpts <- (vec1 + vec2)/2
          widths <- diff(seq.points)
          out1=cumsum(widths*(midpts[-1]))
	out2=sum(widths*(midpts[-1]))
	return(out1/out2)
}


gen.bcond <- function(D, cdf.out, seq.points){
	
	U.cond <- runif(D)	
	b.gen <- apply(U.cond <= cdf.out, 1, function(x){ seq.points[min(which(x))]})
	b.gen[is.na(b.gen)] <- max(seq.points)
	b.gen
}


par.predJW=function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XBhat.update1, Gs, lxN, use.cl= TRUE, Rb = 150 , trunc = FALSE){
  seq.points <- qnorm(tauvec, mean = 0, sd = sqrt(sig2bhat))
  lys <- dat.temp$Y
  rhohat.l <- gpdpar[1]; xi.l <- gpdpar[2]; rhohat.u <- gpdpar[3]; xi.u <- gpdpar[4]
    if(bdist == "Laplace"){
	  clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XBhat.update1", 
	 "comp.term", "gpd.dens", "compfden.b1forintplaplace","lXs","dgpd") , envir = environment())
    	cdf.out=t(parSapply(cl, 1:D, FUN = ianumden.mdp2.claplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, 
	sqrt(sig2bhat), XBhat.update1[smc,]))
    }
    if(bdist == "Normal"){
  clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat", "XBhat.update1", 
	"comp.term", "gpd.dens", "compfden.b1forintp" ,"lXs","dgpd"), envir =   environment())
      cdf.out=t(parSapply(cl, 1:D, FUN = ianumden.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,]))
    } 

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


compfden.b1forintp <- function(b1, Dpick, lyst, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
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
	prod(ll.vec)*dnorm(b1, mean = 0, sd = sdb)
}

bootgensimple2trans <- function(sig2bhat.update, GN, XBhat.update1, CNis, nis, D, areafac.pop, rholxilrhouxiu, lam){
  
  ahats.b <- rnorm(D, mean = 0, sd = sqrt(sig2bhat.update ))
  all.y.b <- XBhat.update1 + as.vector(GN%*%ahats.b)
  lyN.b <- areapredraninterp(all.y.b, rholxilrhouxiu, tauvec)	
   
  q25boot  <-  tapply(hinvtransG1(lyN.b, lam), areafac.pop, quantile, prob = 0.25)
  q50boot <-   tapply(hinvtransG1(lyN.b, lam), areafac.pop, quantile, prob = 0.5)
  q75boot <-  tapply(hinvtransG1(lyN.b, lam), areafac.pop, quantile, prob = 0.75)
  
  lys.b <- lyN.b[smc]
  
  list(lys.b, smc, q25boot = q25boot, q50boot = q50boot, q75boot = q75boot)
  
}


par.bootestsimple2FixSortTrans <- function(dat.tempb, X, areafac.samp, D, Gs, tauvec,  Rb, lxN, b.dist, areafac.pop, smc, use.cl, trunc =FALSE, lam){
  
  initpars <- intparfun(dat.tempb, X, areafac.samp,D, Gs,tauvec,  lxN, use.cl)
  
  rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.tempb, tauvec)
  
  XBbetasiguupdate <- par.updatebetasig2bfunSortFix(tauvec,dat.tempb, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl, Rb )
  
  sig2bhatupdate <- XBbetasiguupdate[[3]]
  betahat <- XBbetasiguupdate[[2]]
  XBhatupdate <- XBbetasiguupdate[[1]]
   
  areapred <- predJWTransFixB0Exp(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl, Rb = 150, trunc, lam ,0)
  
  qhat25b <- areapred[[1]]
  qhat50b <- areapred[[2]]
  qhat75b <- areapred[[3]]
  
  list(qhat25b, qhat50b, qhat75b, sig2bhatupdate, betahat, XBhatupdate, rholxilrhouxiu)
  
}

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





par.bootestsimple2FixSortTransRev <- function(dat.tempb, X, areafac.samp, D, Gs, tauvec,  Rb, lxN, b.dist, areafac.pop, smc, use.cl, trunc =FALSE, lam){
  
  initpars <- intparfun(dat.tempb, X, areafac.samp,D, Gs,tauvec,  lxN, use.cl)
  
  rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.tempb, tauvec)
  
  XBbetasiguupdate <- par.updatebetasig2bfunSortFix(tauvec,dat.tempb, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl, Rb )

  
  sig2bhatupdate <- XBbetasiguupdate[[3]]
  betahat <- XBbetasiguupdate[[2]]
  XBhatupdate <- XBbetasiguupdate[[1]]
 
 rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.tempb, tauvec)

    XBbetasiguupdate <- par.updatebetasig2bfunSortFix(tauvec,dat.tempb, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl, Rb )

  
  sig2bhatupdate <- XBbetasiguupdate[[3]]
  betahat <- XBbetasiguupdate[[2]]
  XBhatupdate <- XBbetasiguupdate[[1]]
 
 rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.tempb, tauvec)
  
  areapred <- predJWTransFixB0Exp(tauvec, dat.tempb, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl, Rb = 150, trunc, lam ,0)
  
  qhat25b <- areapred[[1]]
  qhat50b <- areapred[[2]]
  qhat75b <- areapred[[3]]
  
  list(qhat25b, qhat50b, qhat75b, sig2bhatupdate, betahat, XBhatupdate, rholxilrhouxiu)
  
}





