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



comp.term <- function(index, pred.mat, get.tau.u, lys, tauvec, rhohat.l, xi.l, rhohat.u, xi.u){
	whichtau <- get.tau.u[index]
	if(whichtau == 1){
		lij <- (pred.mat[index,2] + pred.mat[index,1])/2
		dev <- lij - lys[index]
		term <- gpd.dens(dev,rhohat.l, xi.l, tauvec[2], tauvec[1])*(tauvec[1] + tauvec[2])/2		
	}
	if(whichtau > 1 & whichtau != Inf){
		term.num <- tauvec[whichtau] - tauvec[whichtau - 1]
		term.den <- (pred.mat[index,whichtau] - pred.mat[index,whichtau - 1])
		term <- term.num/term.den
	}
	if(whichtau == Inf){
		uij <- (pred.mat[index,length(tauvec)] + pred.mat[index,length(tauvec)-1])/2
		dev <- lys[index] - uij
		term <- (1 - (tauvec[length(tauvec)] + tauvec[length(tauvec)-1])/2)*gpd.dens(dev, rhohat.u, xi.u, tauvec[length(tauvec)-1], tauvec[length(tauvec)])
	}
	term
}


ianumden.mvdpfix=function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
  mat=sapply( seq.points, compfall.b1forintp,  Dpick= d, lys=lys, betahats =betahats.update, tauvec=tauvec, rhohat.l=rhohat.l, xi.l=xi.l, rhohat.u=rhohat.u, xi.u=xi.u, areafac.pop=areafac.pop, smc = smc, sdb=sdb, predfix=predfix)
  mat2=cbind(0, mat[,-dim(mat)[2]])
  midpts <- (mat + mat2)/2
  widths <- diff(seq.points)
  out=sapply(1:3, function(c) sum(widths*midpts[c,-1]))
}


ianumden.mvdpfixlaplace=function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
  mat=sapply( seq.points, compfall.b1forintplaplace,  Dpick= d, lys=lys, betahats =betahats.update, tauvec=tauvec, rhohat.l=rhohat.l, xi.l=xi.l, rhohat.u=rhohat.u, xi.u=xi.u, areafac.pop=areafac.pop, smc = smc, sdb=sdb, predfix=predfix)
  mat2=cbind(0, mat[,-dim(mat)[2]])
  midpts <- (mat + mat2)/2
  widths <- diff(seq.points)
  out=sapply(1:3, function(c) sum(widths*midpts[c,-1]))
}

compfnum.mb1forintp <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- predfix[id,] + b1
	 if(length(id) == 1){
          pred.mat <- matrix(pred.mat, nrow = 1)
        }
        if(length(id) == 1){
          get.tau.u <- min(which(lys[id] <= pred.mat))
        }else{
          get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
        }
	ll.vec <- sapply(1:length(lys[id]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dnorm(b1, mean = 0, sd = sdb)*b1
}

compfall.b1forintplaplace=function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
id <- which(areafac.pop[smc] == Dpick)
pred.mat <- predfix[id,] + b1
get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)

den.b1=prod(ll.vec)*dlp(b1, mean = 0, sd = sdb)
num.mb1=den.b1*b1
num.vb1=num.mb1*b1
return(c(num.mb1, num.vb1, den.b1))
}


compfall.b1forintp=function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
id <- which(areafac.pop[smc] == Dpick)
pred.mat <- predfix[id,] + b1
if(length(id) == 1){
          pred.mat <- matrix(pred.mat, nrow = 1)
          get.tau.u <- min(which(lys[id] <= pred.mat))
        }else{
          get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
        }
ll.vec <- sapply(1:length(lys[id]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)

den.b1=prod(ll.vec)*dnorm(b1, mean = 0, sd = sdb)
num.mb1=den.b1*b1
num.vb1=num.mb1*b1

return(c(num.mb1, num.vb1, den.b1))
}



ianumden.mvdpfix=function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
  mat=sapply( seq.points, compfall.b1forintp,  Dpick= d, lys=lys, betahats =betahats.update, tauvec=tauvec, rhohat.l=rhohat.l, xi.l=xi.l, rhohat.u=rhohat.u, xi.u=xi.u, areafac.pop=areafac.pop, smc = smc, sdb=sdb, predfix=predfix)
  mat2=cbind(0, mat[,-dim(mat)[2]])
  midpts <- (mat + mat2)/2
  widths <- diff(seq.points)
  out=sapply(1:3, function(c) sum(widths*midpts[c,-1]))
}


estconstrained <- function(dat.temp, fml, tauvec){
	rqCur <- rq(fml, data = dat.temp, tau = 0.5, method = "fn")
	coefCur <- coef(rqCur)
	fitCur <- fitted(rqCur)
	iter.taus <- which(tauvec == 0.5)
	coefalls <- coefCur
	repeat{
		iter.taus <- iter.taus + 1
		rqUpdate <- rq(fml, data = dat.temp, tau = tauvec[iter.taus ], method = "fnc", R = lXN, r = lXN%*%coefCur)
		coefCur <- coef(rqUpdate)
		coefalls <- cbind(coefalls, coefCur)
		rqCur <- rqUpdate
		print(iter.taus)
		if(iter.taus == length(tauvec)){break}
	}
	rqCur <- rq(fml, data = dat.temp, tau = 0.5, method = "fn")
	coefCur <- coef(rqCur)
	fitCur <- fitted(rqCur)
	iter.taus <- which(tauvec == 0.5)
	repeat{
		iter.taus <- iter.taus - 1
		rqUpdate <- rq(fml, data = dat.temp, tau = tauvec[iter.taus ], method = "fnc", R = -lXN, r = -lXN%*%coefCur)
		coefCur <- coef(rqUpdate)
		coefalls <- cbind( coefCur, coefalls)
		rqCur <- rqUpdate
		print(iter.taus)
		if(iter.taus == 1){break}
	}
	coefalls
}


par.updatebetasig2bfunConstrFix=function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XB.init, Gs, lxN, use.cl= TRUE, Rb = 150 ){
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
  betahats.update <- estconstrained(dat.temp, as.formula("dev~X"),   tauvec)
  XBhat.update1 <- cbind(1, lxN)%*%betahats.update
  list(XBhat.update1, betahats.update, sig2bhat.update, mub.cond, vebcond)
}


par.predJWFix <-  function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XBhat.update1, Gs, lxN, use.cl= TRUE, Rb = 150 , trunc = FALSE, lam){
  seq.points <- qnorm(tauvec, mean = 0, sd = sqrt(sig2bhat))
  lys <- dat.temp$Y
  rhohat.l <- gpdpar[1]; xi.l <- gpdpar[2]; rhohat.u <- gpdpar[3]; xi.u <- gpdpar[4]
     if(bdist == "Laplace"){
      clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XBhat.update1", 
	 "comp.term", "gpd.dens", "compfall.b1forintplaplace","lXs","dgpd") , envir = environment())
    	out=parSapply(cl, 1:D, FUN = ianumden.mvdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
    if(bdist == "Normal"){
     clusterExport(cl, c("seq.points", "lys",  "tauvec", "rhohat.l", "xi.l", "rhohat.u", "xi.u", "areafac.pop", "smc", "sig2bhat","dlp", "XBhat.update1", 
	"compfnum.mb1forintp", "comp.term", "gpd.dens", "compfall.b1forintp" ,"lXs","dgpd"), envir =   environment())
      out=parSapply(cl, 1:D, FUN = ianumden.mvdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    } 

	num.mub=out[1,]
	num.vb=out[2,]
	den.mub=out[3,]

 
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


betahat.glsfun <- function(sigma2vec, r.ols, k, X.mat, response){
    m <- length(r.ols)
    ci <- 1/sqrt(sigma2vec)/sum(1/sqrt(sigma2vec))
    sig2.tildeu <- sum(ci*(m/(m-k)*r.ols^2 - sigma2vec))
    Vhat1 <- sum(ci^2*((m/(m-k)*r.ols^2 - sigma2vec)-sig2.tildeu)^2)
    sig2u <- max(c(0.5*sqrt(Vhat1), sig2.tildeu))
    W.mat <- diag(1/(sig2u + sigma2vec))
    betahat.gls <- solve(t(X.mat)%*%W.mat%*%X.mat)%*%(t(X.mat)%*%W.mat%*%response)
    list(betahat.gls, sig2u, sig2.tildeu)
}

gpdxi <- function(xi, dvec, rho, taul, tauu){
	prod(sapply(dvec, gpd.dens,  rho, xi, taul, tauu))
}
	
gpd.dens <- function(d, rho, xi, taul, tauu){
	dgpd(d, loc = 0, scale= rho, shape = xi)
}


dlp <- function(t, mean, sd){
	b <- sqrt(sd^2/2)
	exp(-abs(t - mean)/b)/2/b
}