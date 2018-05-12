
predJWTransFixB0Exp <- function(tauvec, dat.temp, sig2bhat, bdist, betahat, gpdpar, areafac.pop, smc, XBhat.update1, Gs, lxN, use.cl= TRUE, Rb = 150 , trunc = FALSE, lam, B = 0){
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
     outall <- sapply( 1:D, FUN =ianumden.mvdpfixlaplace , seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#    num.mubexp <-sapply( 1:D, FUN = ianum.mdplaplacefixexp, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#     num.vb <-sapply( 1:D, FUN = ianum.vdpfixlaplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#    cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
#     den.mub <-  sapply( 1:D,  FUN = iaden.mdp2laplace, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    }
    if(bdist == "Normal"){
      outall <- sapply( 1:D, FUN = ianumden.mdvpfix , seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      num.mubexp  <-sapply( 1:D, FUN =ianum.mdpfixexp, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      num.mub  <-sapply( 1:D, FUN = ianum.mdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#      num.vb <-sapply( 1:D, FUN = ianum.vdpfix, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,])
#     cdf.out1 <- t(sapply( 1:D, FUN = ianum.mdp2.c, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat),XBhat.update1[smc,]))
#      den.mub <-  sapply( 1:D,  FUN = iaden.mdp2, seq.points, lys, betahat, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XBhat.update1[smc,])
    } 
  } 

	 num.mubexp <- outall[4,]
	 num.mub <- outall[1,]
	 num.vb <- outall[2,]
       den.mub <- outall[3,]  

#	print(paste("Len", length(den.mub))

  mub.condexp <- num.mubexp/den.mub
  mub.cond <- num.mub/den.mub
  vebcond <- num.vb/den.mub
    
  all.q.Jtrans <- exp(XBhat.update1)*as.vector(GN%*%mub.condexp)
  all.q.J <-  XBhat.update1 +   as.vector(GN%*%mub.cond)
  t.min.J <- min(all.q.J)
  t.max.J <- max(all.q.J)
  t.seq.J <- seq(t.min.J, t.max.J, length = 100)
  
  if(trunc){
    
    all.q.J[all.q.J <= 0] <- 0
    
  }
  
 # all.q.Jtrans <- t(apply(all.q.J, 1, hinvtransG1, lam) - B)
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




ianumden.mdvpfix=function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
  mat=sapply( seq.points, compfall.b1forintp,  Dpick= d, lys=lys, betahats =betahats.update, tauvec=tauvec, rhohat.l=rhohat.l, xi.l=xi.l, rhohat.u=rhohat.u, xi.u=xi.u, areafac.pop=areafac.pop, smc = smc, sdb=sdb, predfix=predfix)
  mat2=cbind(0, mat[,-dim(mat)[2]])
  midpts <- (mat + mat2)/2
  widths <- diff(seq.points)
  out=sapply(1:4, function(c) sum(widths*midpts[c,-1]))
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
num.mb1exp <- den.b1*exp(b1)
return(c(num.mb1, num.vb1, den.b1, num.mb1exp))
}


ianumden.mvdpfixlaplace=function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
  mat=sapply( seq.points, compfall.b1forintplaplace,  Dpick= d, lys=lys, betahats =betahats.update, tauvec=tauvec, rhohat.l=rhohat.l, xi.l=xi.l, rhohat.u=rhohat.u, xi.u=xi.u, areafac.pop=areafac.pop, smc = smc, sdb=sdb, predfix=predfix)
  mat2=cbind(0, mat[,-dim(mat)[2]])
  midpts <- (mat + mat2)/2
  widths <- diff(seq.points)
  out=sapply(1:4, function(c) sum(widths*midpts[c,-1]))
}


compfall.b1forintplaplace=function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
id <- which(areafac.pop[smc] == Dpick)
pred.mat <- predfix[id,] + b1
get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)

den.b1=prod(ll.vec)*dlp(b1, mean = 0, sd = sdb)
num.mb1=den.b1*b1
num.mb1exp <- den.b1*exp(b1)
num.vb1=num.mb1*b1
return(c(num.mb1, num.vb1, den.b1, num.mb1exp ))
}










