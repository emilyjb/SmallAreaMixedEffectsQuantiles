
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


areapredraninterp  <- function(all.q.Jpop, rholxilrhouxiu, tauvec){
	tausRan <- runif(nrow(all.q.Jpop))
	lowerLimTau <- mean(tauvec[c(1,2)])
	upperLimTau <- mean(tauvec[c(length(tauvec),length(tauvec)-1)])
	tausRanLower <- tausRan[tausRan < lowerLimTau]/lowerLimTau
	tausRanUpper <- (tausRan[tausRan > upperLimTau]  - (tauvec[length(tauvec)]+tauvec[length(tauvec)-1])/2)/(1 - (tauvec[length(tauvec)]+tauvec[length(tauvec)-1])/2)
	ystarlower <- -qgpd(tausRanLower, rholxilrhouxiu[1], rholxilrhouxiu[2]) + apply(all.q.Jpop[which(tausRan < lowerLimTau),c(1,2)], 1, mean)
	ystarupper <-  qgpd(tausRanUpper, rholxilrhouxiu[3], rholxilrhouxiu[4]) + apply(all.q.Jpop[which(tausRan >upperLimTau),c(length(tauvec)-1,length(tauvec))], 1, mean)
	ks <- sapply(tausRan, function(t){ min(which(t <= tauvec))})
	ystar1 <- sapply(1:nrow(all.q.Jpop),interpkreverse,  ks, tausRan, tauvec,all.q.Jpop)
	ystar1[ tausRan < lowerLimTau] <- ystarlower
	ystar1[ tausRan > upperLimTau] <- ystarupper
	ystar1
}


 

areapredraninterp2Old <- function(all.q.Jpop, rholxilrhouxiu, tauvec){
  tausRan <- runif(nrow(areapred$all.q.J))
  lowerLimTau <- mean(tauvec[c(1,2)])
  upperLimTau <- mean(tauvec[c(length(tauvec),length(tauvec)-1)])
  tausRanLower <- tausRan[tausRan < lowerLimTau]/lowerLimTau
  tausRanUpper <- (tausRan[tausRan > upperLimTau]  - (tauvec[length(tauvec)]+tauvec[length(tauvec)-1])/2)/(1 - (tauvec[length(tauvec)]+tauvec[length(tauvec)-1])/2)
  ystarlower <- -qgpd(tausRanLower, rholxilrhouxiu[1], rholxilrhouxiu[2]) + apply(areapred$all.q.J[which(tausRan < lowerLimTau),c(1,2)], 1, mean)
  ystarupper <-  qgpd(tausRanUpper, rholxilrhouxiu[3], rholxilrhouxiu[4]) + apply(areapred$all.q.J[which(tausRan >upperLimTau),c(length(tauvec)-1,length(tauvec))], 1, mean)
  ks <- sapply(tausRan, function(t){ min(which(t <= tauvec))})
  ystar1 <- sapply(1:nrow(all.q.Jpop),interpkreverse,  ks, tauvec, areapred$all.q.J)
  ystar1[ tausRan < lowerLimTau] <- ystarlower
  ystar1[ tausRan > upperLimTau] <- ystarupper
  ystar1
}





