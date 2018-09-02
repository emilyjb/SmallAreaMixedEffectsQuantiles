


####  

htransGJ1 <- function(lys, lam1){
  if(lam1 != 0){
    1/(2*lam1)*(lys^lam1 - 1/lys^lam1)
  }else{
    log(lys)
  }
}


hinvtransG1 <- function(lts, lam1){
  if(lam1 != 0){
	exp(  log( lam1*lts + sqrt( (lts*lam1)^2 + 1))/lam1)
  }else{
	exp(lts)
  }
}

sselam <- function( lam, y, x,  htransGJ1){
	ts <- htransGJ1( y , lam)
	fit <- lm(y~x)$fitted
     	sum( (ts - fit)^2)
}



sselam2 <- function( lambet, y, x,  htransGJ1){
	ts <- htransGJ1( y , lambet[1])
     	sum( abs(ts - lambet[1] - lambet[2]))
}



llhood.max.norm.lam1 <- function(lam1, sigbhat, seq.points, lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  XB.init, htrans){
  lystrans <- htrans(lys, lam1)
  den.mubt <- sapply(1:D, iaden.mdp2, seq.points, lystrans, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sigbhat, XB.init)
  -sum(log(den.mubt))
  
}


llhood.max.laplace.lam1 <- function(lam1, sigbhat, seq.points, lys, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  XB.init, htrans){
  lystrans <- htrans(lys, lam1)
  den.mubt <- sapply(1:D, iaden.mdp2laplace, seq.points, lystrans, coef, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sigbhat, XB.init)
  -sum(log(den.mubt))
  
}







