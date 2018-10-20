##########  Output for initial table of MC properties of predictors:
maxcnt <- 200


outmse1 <- cbind(tapply(apply( mhb25Exps[1:maxcnt,], 2, mean), nis[1:D[1]], mean)*100,
tapply(apply( mhb25s[1:maxcnt,], 2, sd), nis[1:D[1]], mean)/sqrt(nrow(mhb25s[1:maxcnt,])*20)*100,

tapply(apply( uci25Exps[1:maxcnt,] >= qhats.pop25[1:maxcnt,]  & lci25Exps[1:maxcnt,] <= qhats.pop25[1:maxcnt,] , 2, mean), nis[1:D[1]], mean),

tapply(apply((qhats.pop25[1:maxcnt,]  - exp(w25JCs[1:maxcnt,]))^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25[1:maxcnt,]  - EB.norm.25s[1:maxcnt,])^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25[1:maxcnt,]  - q25bs[1:maxcnt,])^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25[1:maxcnt,]  -  qhats.samp25[1:maxcnt,])^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25[1:maxcnt,]  - exp(w25JCs[1:maxcnt,])) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25[1:maxcnt,]  - EB.norm.25s[1:maxcnt,]) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25[1:maxcnt,]  - q25bs[1:maxcnt,]) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25[1:maxcnt,]  -  qhats.samp25[1:maxcnt,]) , 2, mean), nis[1:D[1]], mean)*100 

)



outmse2 <- cbind(tapply(apply( mhb50Exps[1:maxcnt,], 2, mean), nis[1:D[1]], mean)*100,
tapply(apply( mhb50s[1:maxcnt,], 2, sd), nis[1:D[1]], mean)/sqrt(nrow(mhb25s[1:maxcnt,])*20)*100,
tapply(apply( uci50Exps[1:maxcnt,] >= qhats.pop5[1:maxcnt,]  & lci50Exps[1:maxcnt,] <= qhats.pop5[1:maxcnt,] , 2, mean), nis[1:D[1]], mean),
tapply(apply((qhats.pop5[1:maxcnt,]  - exp(w50JCs[1:maxcnt,]))^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5[1:maxcnt,]  - EB.norm.50s[1:maxcnt,])^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5[1:maxcnt,]  - q5bs[1:maxcnt,])^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5[1:maxcnt,]  -  qhats.samp5[1:maxcnt,])^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5[1:maxcnt,]  - exp(w50JCs[1:maxcnt,])) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5[1:maxcnt,]  - EB.norm.50s[1:maxcnt,]) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5[1:maxcnt,] - q5bs[1:maxcnt,]) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5[1:maxcnt,]  - qhats.samp5[1:maxcnt,]) , 2, mean), nis[1:D[1]], mean)*100 
)



 
outmse3 <- cbind(tapply(apply( mhb75Exps[1:maxcnt,], 2, mean), nis[1:D[1]], mean)*100,
tapply(apply( mhb75s[1:maxcnt,], 2, sd), nis[1:D[1]], mean)/sqrt(nrow(mhb25s[1:maxcnt,])*20)*100,
tapply(apply( uci75Exps[1:maxcnt,] >= qhats.pop75[1:maxcnt,]  &   lci75Exps[1:maxcnt,] <=qhats.pop75[1:maxcnt,] , 2, mean), nis[1:D[1]], mean),
tapply(apply((qhats.pop75[1:maxcnt,]  - exp(w75JCs[1:maxcnt,]))^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75[1:maxcnt,]  - EB.norm.75s[1:maxcnt,])^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75[1:maxcnt,]  - q75bs[1:maxcnt,])^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75[1:maxcnt,]  -  qhats.samp75[1:maxcnt,])^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75[1:maxcnt,]  - exp(w75JCs[1:maxcnt,])), 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75[1:maxcnt,]  - EB.norm.75s[1:maxcnt,]), 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75[1:maxcnt,]  - q75bs[1:maxcnt,]), 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75[1:maxcnt,]  -  qhats.samp75[1:maxcnt,]) , 2, mean), nis[1:D[1]], mean)*100 
)


###  Columns of output matrices as follows:
##### MC mean of MSE est.
##### Estimate of MC Standard Error of MC mean of MSE est.  
##### MC coverage
##### MC MSE of LIGPD, NEB, ALD, sample mean
##### MC bias of LIGPD, NEB, ALD, sample mean
rbind(outmse1, outmse2, outmse3)

library("xtable")
xtable(rbind(outmse1, outmse2, outmse3)[,-2], digits = 3)


apply(rbind(outmse1, outmse2, outmse3)[,c(4,1)], 1, diff)/rbind(outmse1, outmse2, outmse3)[,c(4)]*100

RBOutput <- matrix(apply(rbind(outmse1, outmse2, outmse3)[,c(4,1)], 1, diff)/rbind(outmse1, outmse2, outmse3)[,c(4)]*100, ncol = 3, byrow = TRUE)
CovOutput <-  matrix(rbind(outmse1, outmse2, outmse3)[,3], nrow = 3, byrow = TRUE)

cbind(RBOutput, CovOutput)

(rbind(outmse1, outmse2, outmse3)[,8]/100)^2/(rbind(outmse1, outmse2, outmse3)[,4]/100)


