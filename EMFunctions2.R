


#allcdf.d <- sapply(1:D, ianum.mdp.em, seq.points, lys,fit.init$coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat))


dcforem2 <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	sum(sapply(seq.points, compfnum.mb1forintp, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)[-1]*diff(seq.points))
}



ianum.mdp.em <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	cumsum(sapply(seq.points, compfnum.mb1forintp.em, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)[-1]*diff(seq.points))
}



compfnum.mb1forintp.em <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- predfix[id,] + b1
	get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dnorm(b1, mean = 0, sd = sdb)
}

cdf.b.cond <- function(seq.point, D,  lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix, Den){
	num.mub.t <- sapply(1:D, ianum.mdp.em, seq.point, lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)
	num.mub.t/Den
}


gen.bcond <- function(D, cdf.out, seq.points){
	
	U.cond <- runif(D)	
	b.gen <- apply(U.cond <= cdf.out, 1, function(x){ seq.points[min(which(x))]})
	b.gen[is.na(b.gen)] <- max(seq.points)
	b.gen
}

plosscomb <- function(bis, paravec, tauvec, y, X, Gs){
	Btaus <- matrix(paravec, nrow= 2)
	sum(tauvec*sapply(1:length(tauvec),  ploss, tauvec, y,X,Gs,0, bis, Btaus))  
}


ploss1t <-  function(ind,tauvec,y,X,Gs,lambda, bis, Btaus){
    bvec <- c(Btaus, bis)
        tau <- tauvec[ind]
    resid=y-X%*%bvec[1:2] - Gs%*%bvec[-c(1,2)]
    U=tau-as.numeric(resid<0)
    loss=sum(resid*U)
    return(loss)
}


ploss1<- function( paravec,bis, tauvec, y, X, Gs, ind){
	ploss1t(ind, tauvec, y,X,Gs,0, bis, paravec)
}

ploss2 <- function(paravec, bis, tau, y, X, Gs){
    #tau <- tauvec[ind]
    resid=y-X%*%as.vector(paravec) - Gs%*%bis
    U=tau-as.numeric(resid<0)
    loss=sum(resid*U)
    return(loss)
}
	
ploss3 <- function( bis,paravec, tau, y, X, Gs){
    #tau <- tauvec[ind]
    resid=y-X%*%as.vector(paravec) - Gs%*%bis
    U=tau-as.numeric(resid<0)
    loss=sum(resid*U)
    return(loss)
}
	
condmean.ploss1 <- function(paravec, bcond.ran, tauvec, lys, lXs, Gs, ind){ 
	mean(apply(bcond.ran, 2, ploss1, paravec, tauvec, lys, lXs, Gs, ind))
}

condmean.ploss3 <- function(paravec, bcond.ran, tau, lys, lXs, Gs){ 
	mean(apply(bcond.ran, 2, ploss3, paravec, tau, lys, lXs, Gs))
}

opt.par.tau.em <- function(ind, coef.start.all, bcond.ran, tauvec, lys, lXs, Gs){
	coef.start <- coef.start.all[,ind]
	optim(coef.start, condmean.ploss1, bcond.ran=bcond.ran, tauvec=tauvec, lys=lys, lXs=lXs, Gs=Gs, ind)$par
}


opt.par.tau.em3 <- function(tau, coef.start.all, bcond.ran, tauvec, lys, lXs, Gs){
	ind <- which(tauvec == tau)
	coef.start <- coef.start.all[,ind]
	optim(coef.start, condmean.ploss3, bcond.ran=bcond.ran, tau=tau, lys=lys, lXs=lXs, Gs=Gs)$par
}

optimtau <- function(tau, initpar, bis, y, X, Gs, tauvec){
	ind <- which(tauvec == tau)
	optim(initpar[,ind], ploss2, bis=bis, tau=tau, y=y, X=X, Gs=Gs)$par
}
 
condmean.ploss <- function(paravec,bcond.ran, tauvec, lys, lXs, Gs){ 
	mean(apply(bcond.ran, 2, plosscomb, paravec, tauvec, lys, lXs, Gs))
}





