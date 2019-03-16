
tauvec <- seq(0.01, 0.99, by = 0.01)
areafac.samp <- areafac.pop[smc]
Gs <- GN[smc,]
X <- lxN[smc]

######## Obtain initial values of regression parameters as described in Appendix 1
initpars <- intparfunSort(dat.temp, X, areafac.samp,D, Gs,tauvec, lxN, use.cl = FALSE)

######## Store initial values of regression parameters
bj00Cs <- rbind(bj00Cs, initpars$beta[1,])
bj01Cs <- rbind(bj01Cs, initpars$beta[2,])

######## Compute initial values of parameters of extreme value distribution
rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.temp, tauvec)

######## Update estimates of regression parameters
XBbetasiguupdate <- par.updatebetasig2bfunConstrFix(tauvec,dat.temp, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl= FALSE, Rb = 1500 )
sig2bhatupdate <- XBbetasiguupdate[[3]]
betahat <- XBbetasiguupdate[[2]]
XBhatupdate <- XBbetasiguupdate[[1]]

######## Update the estimates of the parameters of the extreme value distribution
rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.temp, tauvec)

######## Update estimates of regression parameters
XBbetasiguupdate <- par.updatebetasig2bfunConstrFix(tauvec,dat.temp, sig2bhatupdate, b.dist, betahat,rholxilrhouxiu , areafac.pop, smc, XBhatupdate[smc,], Gs, lxN, use.cl= FALSE, Rb = 1500 )
sig2bhatupdate <- XBbetasiguupdate[[3]]
betahat <- XBbetasiguupdate[[2]]
XBhatupdate <- XBbetasiguupdate[[1]]

######## Update the estimates of the parameters of the extreme value distribution
rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.temp, tauvec)


######## Store the estimates of the regression parameters 
bj10Cs <- rbind(bj10Cs, betahat[1,])
bj11Cs <- rbind(bj11Cs, betahat[2,])

########  Construct the LIGPD predictors
areapred <- par.predJWFix(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl= FALSE, Rb = 1500 )

#######  LIGPD approximation for the population
all.q.Jpop=areapred$all.q.J

sig2blls <- c(sig2blls, sig2bhatupdate)

bhatCs <- rbind(bhatCs, areapred$mubJ )


#alpha <- -5
#delta <- alpha/sqrt(1+alpha^2)
#w <- sqrt(1/(1-2*delta^2/pi))
#xi <- -w*delta*sqrt(2/pi)


#beta0tvs <- qsn(tauvec, xi, w, alpha)*sig2le + beta0
#beta1tvs <- rep(beta1, length(tauvec))
#XBhattv <- cbind(1, lxN)%*%rbind(beta0tvs, beta1tvs)
#areapredTP <- par.predJWFix(tauvec, dat.temp, sig2lu,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhattv , Gs, lxN,use.cl= FALSE, Rb = 1500 )


 