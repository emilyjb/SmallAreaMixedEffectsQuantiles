


#### Compute the term in the likelihood:

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

compfobs.b <- function(bvec, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc){
	pred.mat <- lXs%*%betahats + kronecker(as.vector(Gs%*%bvec), matrix(rep(1,length(tauvec)),nrow = 1))
	get.tau.u <- apply(lys <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[,1]), comp.term, pred.mat,get.tau.u, lys, tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	tapply(ll.vec, areafac.pop[smc], prod)
}



compfobs.bp <- function(bvec, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, predfix){
	pred.mat <- predfix + kronecker(as.vector(Gs%*%bvec), matrix(rep(1,length(tauvec)),nrow = 1))
	get.tau.u <- apply(lys <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[,1]), comp.term, pred.mat,get.tau.u, lys, tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	tapply(ll.vec, areafac.pop[smc], prod)
}



compfobs.b1p <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, predfix){
	id <- which(areafac.pop[smc == Dpick])
	pred.mat <- predfix[id,] + b1
	get.tau.u <- apply(lys <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys, tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)
}



compfnum.b1 <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- lXs[id,]%*%betahats + b1
	get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dnorm(b1, mean = 0, sd = sdb)
}


compfnum.b1p <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- predfix[id,] + b1
	get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dnorm(b1, mean = 0, sd = sdb)
}


compfnum.mb1forint <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- lXs[id,]%*%betahats + b1
	get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dnorm(b1, mean = 0, sd = sdb)*b1
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



compfnum.vb1forintp <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
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
	prod(ll.vec)*dnorm(b1, mean = 0, sd = sdb)*b1^2
}


dlp <- function(t, mean, sd){
	b <- sqrt(sd^2/2)
	exp(-abs(t - mean)/b)/2/b
}

qlp <- function(p, mean, sd){
	b <- sqrt(sd^2/2)
	if(p <= 0.5){
		b*log(2*p) + mean
	}else{
		mean-b*log(2*(1-p))
	}
}



compfnum.mb1forintplaplace <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- predfix[id,] + b1
	get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dlp(b1, mean = 0, sd = sdb)*b1
}



compfnum.vb1forintplaplace <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- predfix[id,] + b1
	get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lys[id]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dlp(b1, mean = 0, sd = sdb)*b1^2
}



compfnum.mb1forintpsn <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix, snpar){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- predfix[id,] + b1
	get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dsn(b1, xi = snpar[2], omega = snpar[1], alpha = snpar[3])*b1
}



compfnum.vb1forintpsn <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix, snpar){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- predfix[id,] + b1
	get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dsn(b1, xi = snpar[2], omega = snpar[1], alpha = snpar[3])*b1^2
}


#Dpick <- 7
#betahats <- betahats.init
#b1 <- seq.points[1]
#lys <- lys[id]


compfden.b1 <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- lXs[id,]%*%betahats + b1
	get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec) 
}


compfden.b1p <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- predfix[id,] + b1
	get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec) 
}


compfden.b1forint <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- lXs[id,]%*%betahats + b1
	get.tau.u <- apply(lys[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(lXs[id,][,1]), comp.term, pred.mat,get.tau.u, lys[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dnorm(b1, mean = 0, sd = sdb)
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


compfden.b1forintplaplace <- function(b1, Dpick, lyst, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- predfix[id,] + b1
	get.tau.u <- apply(lyst[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(id), comp.term, pred.mat,get.tau.u, lyst[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dlp(b1, mean = 0, sd = sdb)
}


compfden.b1forintpsn <- function(b1, Dpick, lyst, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix,snpar){
	id <- which(areafac.pop[smc] == Dpick)
	pred.mat <- predfix[id,] + b1
	get.tau.u <- apply(lyst[id] <= pred.mat, 1, function(x){ min(which(x))})
	ll.vec <- sapply(1:length(id), comp.term, pred.mat,get.tau.u, lyst[id], tauvec,  rhohat.l, xi.l, rhohat.u, xi.u)
	prod(ll.vec)*dsn(b1, xi = snpar[2], omega = snpar[1], alpha = snpar[3])
}




ianum.md <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
	sum(sapply(seq.points, compfnum.mb1forint, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb)[-1]*diff(seq.points))
}

iaden.md <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
	sum(sapply(seq.points, compfden.b1forint, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb)[-1]*diff(seq.points))
}


ianum.mdp <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
  vec1 <- sum(sapply( seq.points, compfnum.mb1forintp,  Dpick= d, lys=lys, betahats =betahats.update, tauvec=tauvec, rhohat.l=rhohat.l, xi.l=xi.l, rhohat.u=rhohat.u, xi.u=xi.u, areafac.pop=areafac.pop, smc = smc, sdb=sdb, predfix=predfix))
  vec2 <- c(0, vec1[-length(vec1)])
  midpts <- (vec1 + vec2)/2
  widths <- diff(seq.points)
  cumsum(widths*(midpts[-1]))
	#sum(sapply(seq.points, compfnum.mb1forintp, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)[-1]*diff(seq.points))
}


ianum.vdp <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
  vec1 <-  sum(sapply(seq.points, compfnum.vb1forintp, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)[-1]*diff(seq.points))
  vec2 <- c(0, vec1[-length(vec1)])
  midpts <- (vec1 + vec2)/2
  widths <- diff(seq.points)
  cumsum(widths*(midpts[-1]))
}



iaden.mdp <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	sum(sapply(seq.points, compfden.b1forintp, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)[-1]*diff(seq.points))
}





ianum.mdplaplace <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	sum(sapply(seq.points, compfnum.mb1forintplaplace, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)[-1]*diff(seq.points))
}


ianum.vdplaplace <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	sum(sapply(seq.points, compfnum.vb1forintplaplace, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)[-1]*diff(seq.points))
}



iaden.mdplaplace <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix){
	sum(sapply(seq.points, compfden.b1forintplaplace, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix)[-1]*diff(seq.points))
}


ianum.mdpsn <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix, snpar){
	sum(sapply(seq.points, compfnum.mb1forintpsn, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix, snpar)[-1]*diff(seq.points))
}


ianum.vdpsn <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix,snpar){
	sum(sapply(seq.points, compfnum.vb1forintpsn, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix, snpar)[-1]*diff(seq.points))
}



iaden.mdpsn <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix,snpar){
	sum(sapply(seq.points, compfden.b1forintpsn, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb, predfix, snpar)[-1]*diff(seq.points))
}



#####################################   Write function to smooth the XB-hat-init     #####################################

smoothfit.row0 <- function(row, XBhatfix){
	testdif <- t(apply(XBhatfix, 1, diff))
	index.0 <- which(testdif[row,] == 0)
	if( 1 %in% index.0){
		XBhatfix[row,1] <- XBhatfix[row,1] - testdif[min(which(testdif>0))]
		testdif <- t(apply(XBhatfix, 1, diff))
		index.0 <- which(testdif[row,] == 0)
	}
	zc <- length(index.0) > 0
	yo <- XBhatfix[row,]
	while(zc){			
		ynew <- yo
		ynew[index.0] <- yo[index.0-1] + (tauvec[index.0] - tauvec[index.0-1])*(yo[index.0+1] - yo[index.0-1])/(tauvec[index.0+1] - tauvec[index.0-1])
		yo <- ynew
		testdif <- diff(yo)
		index.0 <- which(testdif == 0)
		zc <- length(index.0) > 0
		#ynew <- yo
		#ynew[index.0] <- yo[index.0-1] + (tauvec[index.0] - tauvec[index.0-1])*(yo[index.0+1] - yo[index.0-1])/(tauvec[index.0+1] - tauvec[index.0-1])
		#diff(ynew)
	}
	yo
}


smoothfit.row <- function(row, XBhatfix, tauvec){
testdif <- t(apply(XBhatfix, 1, diff))
index.0 <- which(testdif[row,] == 0)
if(length(index.0) > 0){
if( 1 %in% index.0){
XBhatfix[row,1] <- XBhatfix[row,1] - testdif[row,][min(which(testdif[row,]>0))]
testdif <- t(apply(XBhatfix, 1, diff))
index.0 <- which(testdif[row,] == 0)
}
if( length(tauvec) %in% index.0){
XBhatfix[row,length(tauvec)] <- XBhatfix[row,length(tauvec)] +  testdif[row,][max(which(testdif[row,]>0))]
testdif <- t(apply(XBhatfix, 1, diff))
index.0 <- which(testdif[row,] == 0)
}
group <- rep(0,length(index.0))
group[1] <- 1
a.gi <- 0
if(length(index.0) > 1){
di0 <- diff(index.0)
repeat{
a.gi <- a.gi + 1
if(di0[a.gi] == 1){
group[a.gi+1] <- group[a.gi]
}else{
group[a.gi+1] <- group[a.gi] + 1 
}
if(a.gi == length(di0)){break}
}
}
yo <- XBhatfix[row,]
iter.g <- 0
repeat{
iter.g <- iter.g + 1
whichind <- index.0[which(group == iter.g)]
yvec <- yo[whichind]
yl <- yo[min(whichind)-1]
yu <- yo[max(whichind)+1]
yo[whichind] <- yl + (tauvec[whichind] - tauvec[min(whichind)-1])*(yu - yl)/(tauvec[max(whichind)+1] - tauvec[min(whichind)-1])
if(iter.g == max(group)){break}
}
yo
}else{
XBhatfix[row,]
}
}


##### XBNsm <- sapply(1:nrow(XBhatfix), smoothfit.row, XBhatfix)	
##### iter.sm <- 0
##### repeat{
#####	iter.sm <- iter.sm + 1
#####	smout <- smoothfit.row(iter.sm, XBhatfix)
#####	print(paste(iter.sm))
#####}



#isoint <- isoreg(tauvec, fit.init$coef[1,])

#lines(fit.init$coef[1,])

#compfcond.b1 <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
#	num <- compfnum.b1(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  sdb)
#	den <- compfden.b1(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  sdb)
#	num/den
#}



#compfcond.b1forint <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
#	num <- compfnum.b1(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  sdb)
#	den <- compfden.b1forint(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  sdb)
#	num/den
#}#


#sapply(seq.points, compfnum.mb1forint, 7,  lys, betahats.init, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2bhat))


#compfcond.mb1 <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
#	num <- compfnum.b1(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  sdb)
#	den <- compfden.b1(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  sdb)
#	num*b1/den
#}


#compfcond.mb1forint <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
#	num <- compfnum.b1(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  sdb)
#	den <- compfden.b1forint(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  sdb)
#	num*b1/den
#}



#iacond.d <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
#	cumsum(sapply(seq.points, compfcond.b1forint, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb)[-1]*diff(seq.points))
#}


#iamean.d <- function(d, seq.points, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
#	sum(sapply(seq.points, compfcond.mb1forint, d, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb)[-1]*diff(seq.points))
#}#


#compfcond.mb1 <- function(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sdb){
#	num <- compfnum.b1(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  sdb)
#	den <- compfden.b1(b1, Dpick, lys, betahats, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc,  sdb)
#	num*b1/den
#}




#compfcond.mb1forint(0.5, 1, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2uhatmc))
#compfcond.mb1forint(0.5, 2, lys, betahats.update, tauvec, rhohat.l, xi.l, rhohat.u, xi.u, areafac.pop, smc, sqrt(sig2uhatmc))

snd <- function(del, V, bs){
	omega <- sqrt(V/(1 -2*del^2/pi))
	xi <- -omega*del*sqrt(2/pi)
	alp <- sqrt(del^2/(1-del^2))
	-sum(sapply(bs, dsn, xi=xi, omega=omega, alpha = alp, log = TRUE) )
}


parfun.del <- function(del, V, bs){
	omega <- sqrt(V/(1 -2*del^2/pi))
	xi <- -omega*del*sqrt(2/pi)
	alp <- sqrt(del^2/(1-del^2))
	c(omega, xi, alp)
}

	


	




