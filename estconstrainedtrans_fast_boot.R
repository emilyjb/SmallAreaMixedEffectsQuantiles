
tauvec <- seq(0.01, 0.99, by = 0.01)
areafac.samp <- areafac.pop[smc]
Gs <- GN[smc,]
X <- lxN[smc]

lamcur <- 0
yobs <- hinvtransG1(lys, lamcur)
yst <- htransGJ1(yobs, lamcur)
 
dat.temp <- data.frame(Y = yst, X) 

initpars <- intparfunSort(dat.temp, X, areafac.samp,D, Gs,tauvec, lxN, use.cl = FALSE)
rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.temp, tauvec)

lamcur <- 0
yst <- htransGJ1(yobs, lamcur)
dat.temp <- data.frame(Y = yst, X = X, areafac.pop = areafac.samp)
lamcursels <- c(lamcursels, lamcur)

bj00Cs <- rbind(bj00Cs, initpars$beta[1,])
bj01Cs <- rbind(bj01Cs, initpars$beta[2,])


XBbetasiguupdate <- par.updatebetasig2bfunConstrFix(tauvec,dat.temp, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl= FALSE, Rb = 1500 )

sig2bhatupdate <- XBbetasiguupdate[[3]]
betahat <- XBbetasiguupdate[[2]]
XBhatupdate <- XBbetasiguupdate[[1]]

rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.temp, tauvec)


bj10Cs <- rbind(bj10Cs, betahat[1,])
bj11Cs <- rbind(bj11Cs, betahat[2,])

sig2blls <- c(sig2blls, sig2bhatupdate)


areapred <- predJWTransFixB0Exp(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl= FALSE, Rb = 1500, trunc = FALSE, lamcur )







