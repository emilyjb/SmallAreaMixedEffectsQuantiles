

### Functions for optimization for random bi
compfden.b1forintp2 <- function(b1, Dpick, lyst, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        id <- which(areafac.pop[smc] == Dpick)
        pred.mat <- predfix[id,] + b1
        get.tau.u <- apply(lyst[id] <= pred.mat, 1, function(x){ min(which(x))})
        ll.vec <- sapply(1:length(id), comp.term, pred.mat,get.tau.u, lyst[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
        prod(ll.vec) 
}


iaden.mdpran <- function(d, L, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	  seq.points <- sort(rnorm(L, mean = 0, sd = sdb))
        sum(sapply(seq.points, compfden.b1forintp2, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)[-1]*diff(seq.points))
}



llhood.max.norm.ran <- function(sig2bhat, L, lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  XB.init){
	
	den.mubt <- sapply(1:D, iaden.mdpran, L, lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init)
	-sum(log(den.mubt))

}



### Functions for optimization with sequence of points 


iaden.mdp2 <-  function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        vec1 <- sapply(seq.points, compfden.b1forintp, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix) 
	  vec2 <- c(0, vec1[-length(vec1)])
	  midpts <- (vec1 + vec2)/2
	  widths <- diff(seq.points)
	  sum(widths*(midpts[-1]))
	
}


iaden.mdp2sn <-  function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix, snpar){
        vec1 <- sapply(seq.points, compfden.b1forintpsn.ll, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix,snpar) 
	  vec2 <- c(0, vec1[-length(vec1)])
	  midpts <- (vec1 + vec2)/2
	  widths <- diff(seq.points)
	  sum(widths*(midpts[-1]))
	
}

 compfden.b1forintpsn.ll <- function(b1, Dpick, lyst, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix,sigdel){
	  V <- sigdel[1]
	  del <- sigdel[2]
        omega <- sqrt(V/(1 -2*del^2/pi))
        xi <- -omega*del*sqrt(2/pi)
        alp <- sqrt(del^2/(1-del^2))
        id <- which(areafac.pop[smc] == Dpick)
        pred.mat <- predfix[id,] + b1
        get.tau.u <- apply(lyst[id] <= pred.mat, 1, function(x){ min(which(x))})
        ll.vec <- sapply(1:length(id), comp.term, pred.mat,get.tau.u, lyst[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
        prod(ll.vec)*dsn(b1, xi = xi, omega = omega, alpha = alp)
}


llhood.max.sn <- function(sigdel, seq.points, lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  XB.init){
	
	den.mubt <- sapply(1:D, iaden.mdp2sn, seq.points, lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sigbhat, XB.init, sigdel)
	-sum(log(den.mubt))

}

#llhood.max.sn(c(sig2bhat, delhat),seq.points, lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  XB.init) 
#iaden.mdp2sn(1, seq.points, lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat), XB.init, c(sig2bhat, delhat))



llhood.max.norm <- function(sigbhat, seq.points, lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  XB.init){
	
	den.mubt <- sapply(1:D, iaden.mdp2, seq.points, lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sigbhat, XB.init)
	-sum(log(den.mubt))

}

#snpar.opt <- optim(sig2bhat, llhood.max.norm.ran,L = 100, lys=lys, coef = fit.init$coef, tauvec = tauvec, rhohat.l = rhohat.l, xi.l = xi.l, rhohat.u =rhohat.u, xi.u=xi.u, areafac.pop =areafac.pop,smc= smc,  XB.init=XBhat.update1)




ianum.mdp2.c <-  function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        vec1 <- sapply(seq.points, compfden.b1forintp, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)
	  vec2 <- c(0, vec1[-length(vec1)])
	  midpts <- (vec1 + vec2)/2
	  widths <- diff(seq.points)
	  cumsum(widths*(midpts[-1]))
	
}


ianum.mdp2.claplace <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        vec1 <- sapply(seq.points, compfden.b1forintplaplace, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)
          vec2 <- c(0, vec1[-length(vec1)])
          midpts <- (vec1 + vec2)/2
          widths <- diff(seq.points)
          cumsum(widths*(midpts[-1]))
        
}

iaden.mdp2laplace <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
        vec1 <- sapply(seq.points, compfden.b1forintplaplace, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix) 
          vec2 <- c(0, vec1[-length(vec1)])
          midpts <- (vec1 + vec2)/2
          widths <- diff(seq.points)
          sum(widths*(midpts[-1]))
        
}


snloss <- function(sigdel, bs){
	  V <- sigdel[1]
	  del <- sigdel[2]
        omega <- sqrt(V/(1 -2*del^2/pi))
        xi <- -omega*del*sqrt(2/pi)
        alp <- sqrt(del^2/(1-del^2))
 	  -sum(sapply(bs, dsn, xi = xi, omega = omega, alpha = alp, log = TRUE))
}

condmeansnloss <- function(sigdel, bis){
	mean(apply(bis, 2, snloss, sigdel))
}


#optim(c(sig2bhat, delhat), condmeansnloss, bis = bcond.ran)

#snloss(bcond.ran[,1], c(sig2bhat, delhat))
#snloss(bcond.ran[,2], c(sig2bhat, delhat))




ianum.mdp2.csn <-  function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix, snpar){
        vec1 <- sapply(seq.points, compfden.b1forintpsn, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix, snpar)
	  vec2 <- c(0, vec1[-length(vec1)])
	  midpts <- (vec1 + vec2)/2
	  widths <- diff(seq.points)
	  cumsum(widths*(midpts[-1]))
	
}


ianum.mdpst <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  predfix, stpar){
        sum(sapply(seq.points, compfnum.mb1forintpst, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  predfix, stpar)[-1]*diff(seq.points))
}




compfnum.mb1forintpst <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  predfix, stpar){
        id <- which(areafac.pop[smc] == Dpick)
        pred.mat <- predfix[id,] + b1
        get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
        ll.vec <- sapply(1:length(lys[id]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
        prod(ll.vec)*dt.scaled(b1, mean = 0, df = stpar[2], sd = stpar[1])*b1
}

compfden.b1forintpst <- function(b1, Dpick, lyst, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  predfix, stpar){
        id <- which(areafac.pop[smc] == Dpick)
        pred.mat <- predfix[id,] + b1
        get.tau.u <- apply(lyst[id] <= pred.mat, 1, function(x){ min(which(x))})
        ll.vec <- sapply(1:length(id), comp.term, pred.mat,get.tau.u, lyst[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
        prod(ll.vec)*dt.scaled(b1, mean = 0, df = stpar[2], sd = stpar[1])
}


ldtscaledprofdf <- function(d,  b, sig){
	-sum(log(dt.scaled(b, mean = 0, df = d, sd = sqrt(sig^2/(d/(d-2))))))
}


ldtscaledprofsig <- function(sig,  b, d){
	-sum(log(dt.scaled(b, mean = 0, df = d, sd = sqrt(sig^2/(d/(d-2))))))
}

ldtscaledprofdfbran <- function(d,  b, sig){
	mean(apply(b, 2, function(x){-sum(log(dt.scaled(x, mean = 0, df = d, sd = sqrt(sig^2/(d/(d-2)))))) }))
}


ldtscaledprofsigbran <- function(sig,  b, d){
	mean(apply(b, 2, function(x){-sum(log(dt.scaled(b, mean = 0, df = d, sd = sqrt(sig^2/(d/(d-2)))))) }))
}


iaden.mdpst <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, predfix, stpar){
        sum(sapply(seq.points, compfden.b1forintpst, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, predfix, stpar)[-1]*diff(seq.points))
}


ianum.mdp2.cst <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  predfix, stpar){
        vec1 <- sapply(seq.points, compfden.b1forintpst, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,   predfix, stpar)
        vec2 <- c(0, vec1[-length(vec1)])
        midpts <- (vec1 + vec2)/2
        widths <- diff(seq.points)
        cumsum(widths*(midpts[-1]))
        
}

iaden.mdp2st <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, predfix, stpar){
        vec1 <- sapply(seq.points, compfden.b1forintpst, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  predfix, stpar) 
        vec2 <- c(0, vec1[-length(vec1)])
        midpts <- (vec1 + vec2)/2
        widths <- diff(seq.points)
        sum(widths*(midpts[-1]))
        
}

iter.opt <- function(bcond.ran, sigd){
	iter.opt <- 0
	par.cur <- sigd
	repeat{
	iter.opt <- iter.opt + 1
	d.opt <- optimize(ldtscaledprofdfbran, c(2.01, 15), bcond.ran, par.cur[1])$minimum
	s.opt <- mean(bcond.ran^2)*D/(D-2)
	par.new <- c(s.opt, d.opt)
	if(max(abs(par.new - par.cur)) < 0.05){break}
	par.cur <- par.new
	if(iter.opt > 5){break}
	}
	par.cur
}


ianum.vdpst <-  function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  predfix, stpar){
        sum(sapply(seq.points, compfnum.vb1forintpst, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, predfix, stpar)[-1]*diff(seq.points))
}

compfnum.vb1forintpst <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, predfix, stpar){
        id <- which(areafac.pop[smc] == Dpick)
        pred.mat <- predfix[id,] + b1
        get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
        ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
        prod(ll.vec)*dt.scaled(b1, mean = 0, df = stpar[2], sd = stpar[1])*b1^2
}

ldtscaledprofdf <- function(d,  b, sig){
	  b <- b/sqrt(sig)
        -sum(log(dt.scaled(b, mean = 0, df = d, sd = 1)))
}





ianum.mdpfix <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
  vec1 <- sapply( seq.points, compfnum.mb1forintp,  Dpick= d, lys=lys, betahats =betahats.update, tauvec=tauvec, rhohat.l=rhohat.l, xi.l=xi.l, rhohat.u=rhohat.u, xi.u=xi.u, areafac.pop=areafac.pop, smc = smc, sdb=sdb, predfix=predfix)
  vec2 <- c(0, vec1[-length(vec1)])
  midpts <- (vec1 + vec2)/2
  widths <- diff(seq.points)
  sum(widths*(midpts[-1]))
  # sum(sapply(seq.points, compfnum.mb1forintp, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)[-1]*diff(seq.points))
}

ianum.vdpfix <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
  vec1 <-  sapply(seq.points, compfnum.vb1forintp, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)
  vec2 <- c(0, vec1[-length(vec1)])
  midpts <- (vec1 + vec2)/2
  widths <- diff(seq.points)
  sum(widths*(midpts[-1]))
}


ianum.mdpfixlaplace <-  function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
   vec1 <-  sapply(seq.points, compfnum.mb1forintplaplace, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)
  vec2 <- c(0, vec1[-length(vec1)])
  midpts <- (vec1 + vec2)/2
  widths <- diff(seq.points)
  sum(widths*(midpts[-1]))
}

ianum.vdpfixlaplace <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
  vec1 <-  sapply(seq.points, compfnum.vb1forintplaplace, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)
  vec2 <- c(0, vec1[-length(vec1)])
  midpts <- (vec1 + vec2)/2
  widths <- diff(seq.points)
  sum(widths*(midpts[-1]))
}


