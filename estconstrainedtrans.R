
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


#lamseq <- seq(0, 0.9, by = 0.1)
#seq.points <- qnorm(tauvec, mean = 0, sd = initpars$sig2bhat)
#lamcurs <- sapply(lamseq,  llhood.max.norm.lam1,  sigbhat = sqrt(initpars$sig2bhat), seq.points = seq.points, lys =  yobs,coef = initpars$beta, tauvec = tauvec, rhohat.l = rholxilrhouxiu[1], xi.l =rholxilrhouxiu[2], rhohat.u = rholxilrhouxiu[3],xi.u= rholxilrhouxiu[4],areafac.pop= areafac.pop, smc = smc, XB.init = initpars$XBinit, htrans = htransGJ1)
lamcur <- 0
yst <- htransGJ1(yobs, lamcur)
dat.temp <- data.frame(Y = yst, X = X, areafac.pop = areafac.samp)
lamcursels <- c(lamcursels, lamcur)

bj00Cs <- rbind(bj00Cs, initpars$beta[1,])
bj01Cs <- rbind(bj01Cs, initpars$beta[2,])


XBbetasiguupdate <- updatebetasig2bfunConstrFix(tauvec,dat.temp, initpars$sig2bhat, b.dist, initpars$beta,rholxilrhouxiu , areafac.pop, smc, initpars$XBinit[smc,], Gs, lxN, use.cl= FALSE, Rb = 1500 )

sig2bhatupdate <- XBbetasiguupdate[[3]]
betahat <- XBbetasiguupdate[[2]]
XBhatupdate <- XBbetasiguupdate[[1]]

rholxilrhouxiu <- gpdparMatch(XBhatupdate[smc,], dat.temp, tauvec)


bj10Cs <- rbind(bj10Cs, betahat[1,])
bj11Cs <- rbind(bj11Cs, betahat[2,])

sig2blls <- c(sig2blls, sig2bhatupdate)


areapred <- predJWTransFixB0Exp(tauvec, dat.temp, sig2bhatupdate,  b.dist, betahat, rholxilrhouxiu , areafac.pop, smc, XBhatupdate, Gs, lxN,use.cl= FALSE, Rb = 1500, trunc = FALSE, lamcur )


w25JCs <- rbind(w25JCs, areapred[[1]])
w50JCs <- rbind(w50JCs, areapred[[2]])
w75JCs <- rbind(w75JCs, areapred[[3]])
w90JCs <- rbind(w90JCs, areapred[[4]])
w10JCs <- rbind(w10JCs, areapred[[5]])


#ystarMC <-replicate(500, areapredraninterp(areapred$all.q.J, rholxilrhouxiu, tauvec))
#ystarMC <- apply(ystarMC, 2,   hinvtransG1, lam)
#ystarMC[smc,] <- lys

#w10JRs <- rbind(w25JRs, sapply(1:D, function(i){ quantile(as.vector( ystarMC[areafac.pop == i,]), prob = 0.1)}))
#w25JRs <- rbind(w25JRs, sapply(1:D, function(i){ quantile(as.vector( ystarMC[areafac.pop == i,]), prob = 0.25)}))
#w50JRs <- rbind(w50JRs, sapply(1:D, function(i){ quantile(as.vector( ystarMC[areafac.pop == i,]), prob = 0.5)}))
#w75JRs <- rbind(w75JRs, sapply(1:D, function(i){ quantile(as.vector( ystarMC[areafac.pop == i,]), prob = 0.75)}))
#w90JRs <- rbind(w90JRs, sapply(1:D, function(i){ quantile(as.vector( ystarMC[areafac.pop == i,]), prob = 0.9)}))



