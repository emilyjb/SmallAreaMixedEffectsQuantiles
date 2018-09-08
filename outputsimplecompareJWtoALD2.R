rm(list = ls(all = TRUE))
#### Load in appropriate Rdata file; for example:

load("SNOutputInterpF1LaplaceFix_temp90.Rdata")

#### Set number of iterations to use (typically, maxcnt <- cnt; useful if terminate simulation before completion)
maxcnt <- (cnt-1)


MSE25=apply((w25JCs[1:maxcnt,] - qhats.pop25[1:maxcnt,])^2, 2, mean)
MSE50=apply((w50JCs[1:maxcnt,] - qhats.pop5[1:maxcnt,])^2, 2, mean)
MSE75=apply((w75JCs[1:maxcnt,] - qhats.pop75[1:maxcnt,])^2, 2, mean)

B.MSE25=colMeans(mhb25s);B.MSE50=colMeans(mhb50s);B.MSE75=colMeans(mhb75s);
RB=rbind(100*((B.MSE25-MSE25)/MSE25),100*((B.MSE50-MSE50)/MSE50),
100*((B.MSE75-MSE75)/MSE75))


CR.f=function(w25JCs, qhats.pop25, mhb25s){
est.q=w25JCs
tr.q=qhats.pop25
b.MSE=mhb25s
L.CI=est.q-1.96*sqrt(b.MSE)
U.CI=est.q+1.96*sqrt(b.MSE)
CR=apply(L.CI <=tr.q & tr.q <=U.CI, 2, mean)
return(CR)
}

CR=rbind(CR.f(w25JCs, qhats.pop25, mhb25s),CR.f(w50JCs, qhats.pop5, mhb50s),
CR.f(w75JCs, qhats.pop75, mhb75s))

RB.Table1=sapply(1:3, function(i) tapply(RB[i,],nis,mean))
CR.Table1=sapply(1:3, function(i) tapply(CR[i,],nis,mean))


