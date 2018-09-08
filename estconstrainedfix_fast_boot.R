
tauvec <- seq(0.01, 0.99, by = 0.01)
areafac.samp <- areafac.pop[smc]
Gs <- GN[smc,]
X <- lxN[smc]

initpars <- intparfun(dat.temp, X, areafac.samp,D, Gs,tauvec, lxN, use.cl = FALSE)

bj00Cs <- rbind(bj00Cs, initpars$beta[1,])
bj01Cs <- rbind(bj01Cs, initpars$beta[2,])

rholxilrhouxiu <- gpdparMatch(initpars$XBinit[smc,], dat.temp, tauvec)
XBbetasiguupdate <- par.updatebetasig2bfunConstrFix(tauvec,dat.temp, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl= FALSE, Rb = 1500 )

sig2bhatupdate <- XBbetasiguupdate[[3]]
betahat <- XBbetasiguupdate[[2]]
XBhatupdate <- XBbetasiguupdate[[1]]
 
bj10Cs <- rbind(bj10Cs, betahat[1,])
bj11Cs <- rbind(bj11Cs, betahat[2,])
areapred <- par.predJWFix(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl= FALSE, Rb = 1500 )

 
all.q.Jpop=areapred$all.q.J

 