
tauvec <- seq(0.01, 0.99, by = 0.01)
areafac.samp <- areafac.pop[smc]
Gs <- GN[smc,]
X <- lxN[smc]

####  Transform data for estimation
lamcur <- 0
yobs <- hinvtransG1(lys, lamcur)
yst <- htransGJ1(yobs, lamcur)
 
dat.temp <- data.frame(Y = yst, X) 

######## Obtain initial values of regression parameters as described in Appendix 1
initpars <- intparfunSort(dat.temp, X, areafac.samp,D, Gs,tauvec, lxN, use.cl = FALSE)

######## Compute initial values of parameters of extreme value distribution
rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.temp, tauvec)

#######  Transformation material repeated here due to previous code where transformation parameter was estimated
lamcur <- 0
yst <- htransGJ1(yobs, lamcur)
dat.temp <- data.frame(Y = yst, X = X, areafac.pop = areafac.samp)
lamcursels <- c(lamcursels, lamcur)

######## Store initial values of regression parameters
bj00Cs <- rbind(bj00Cs, initpars$beta[1,])
bj01Cs <- rbind(bj01Cs, initpars$beta[2,])

######## Update estimates of regression parameters
XBbetasiguupdate <- par.updatebetasig2bfunConstrFix(tauvec,dat.temp, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl= FALSE, Rb = 1500 )

sig2bhatupdate <- XBbetasiguupdate[[3]]
betahat <- XBbetasiguupdate[[2]]
XBhatupdate <- XBbetasiguupdate[[1]]

######## Update the estimates of the parameters of the extreme value distribution
rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.temp, tauvec)



######## Update estimates of regression parameters
XBbetasiguupdate <- par.updatebetasig2bfunConstrFix(tauvec,dat.temp, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl= FALSE, Rb = 1500 )

########  Construct the LIGPD predictors
areapred1 <- predJWTransFixB0Exp(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl= FALSE, Rb = 1500, trunc = FALSE, lamcur )


sig2bhatupdate <- XBbetasiguupdate[[3]]
betahat <- XBbetasiguupdate[[2]]
XBhatupdate <- XBbetasiguupdate[[1]]

######## Update the estimates of the parameters of the extreme value distribution
rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.temp, tauvec)

 

########  Construct the LIGPD predictors
areapred2 <- predJWTransFixB0Exp(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl= FALSE, Rb = 1500, trunc = FALSE, lamcur )


######## Store the estimates of the regression parameters  and variance
bj10Cs <- rbind(bj10Cs, betahat[1,])
bj11Cs <- rbind(bj11Cs, betahat[2,])


sig2blls <- c(sig2blls, sig2bhatupdate)





