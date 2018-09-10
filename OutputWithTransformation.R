
rm(list = ls(all = TRUE))
#### Load in appropriate Rdata file; for example:

load("SNTransBoot3Normal.Rdata")

#### Set maximum number of iterations to use (typically, maxcnt <- cnt)
maxcnt <- cnt 

##########  Average MSE by area size

MCMSETab10 <-rbind( tapply(apply((w10JCs[1:maxcnt,] - qhats.pop10[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((q10bs[1:maxcnt,] - qhats.pop10[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((EB.norm.10s[1:maxcnt,] - qhats.pop10[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((exp(qhats.samp10[1:maxcnt,]) - qhats.pop10[1:maxcnt,])^2, 2, mean), nis, mean))
rownames(MCMSETab10) <- c("LIGPD","ALD","NEB","Direct")

MCMSETab25 <- rbind(tapply(apply((w25JCs[1:maxcnt,] - qhats.pop25[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((q25bs[1:maxcnt,] - qhats.pop25[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((EB.norm.25s[1:maxcnt,] - qhats.pop25[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((exp(qhats.samp25[1:maxcnt,]) - qhats.pop25[1:maxcnt,])^2, 2, mean), nis, mean))
rownames(MCMSETab25) <- c("LIGPD","ALD","NEB","Direct")

MCMSETab50 <- rbind(tapply(apply((w50JCs[1:maxcnt,] - qhats.pop5[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((q5bs[1:maxcnt,] - qhats.pop5[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((EB.norm.50s[1:maxcnt,] - qhats.pop5[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((exp(qhats.samp5[1:maxcnt,]) - qhats.pop5[1:maxcnt,])^2, 2, mean), nis, mean))
rownames(MCMSETab50) <- c("LIGPD","ALD","NEB","Direct")

MCMSETab75 <- rbind(tapply(apply((w75JCs[1:maxcnt,] - qhats.pop75[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((q75bs[1:maxcnt,] - qhats.pop75[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((EB.norm.75s[1:maxcnt,] - qhats.pop75[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((exp(qhats.samp75[1:maxcnt,]) - qhats.pop75[1:maxcnt,])^2, 2, mean), nis, mean))
rownames(MCMSETab75) <- c("LIGPD","ALD","NEB","Direct")

MCMSETab90 <- rbind(tapply(apply((w90JCs[1:maxcnt,] - qhats.pop90[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((q90bs[1:maxcnt,] - qhats.pop90[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((EB.norm.90s[1:maxcnt,] - qhats.pop90[1:maxcnt,])^2, 2, mean), nis, mean),
	tapply(apply((exp(qhats.samp90[1:maxcnt,]) - qhats.pop90[1:maxcnt,])^2, 2, mean), nis, mean))
rownames(MCMSETab90) <- c("LIGPD","ALD","NEB","Direct")



##########  Average Bias by area size

MCBiasTab10 <- rbind(tapply(apply((w10JCs[1:maxcnt,] - qhats.pop10[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((q10bs[1:maxcnt,] - qhats.pop10[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((EB.norm.10s[1:maxcnt,] - qhats.pop10[1:maxcnt,]), 2, mean), nis, mean),
	tapply(apply((exp(qhats.samp10[1:maxcnt,]) - qhats.pop10[1:maxcnt,])^2, 2, mean), nis, mean))
rownames(MCBiasTab10) <- c("LIGPD","ALD","NEB","Direct")

MCBiasTab25 <- rbind(tapply(apply((w25JCs[1:maxcnt,] - qhats.pop25[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((q25bs[1:maxcnt,] - qhats.pop25[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((EB.norm.25s[1:maxcnt,] - qhats.pop25[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((exp(qhats.samp25[1:maxcnt,] ) - qhats.pop25[1:maxcnt,]), 2, mean), nis, mean))
rownames(MCBiasTab25) <- c("LIGPD","ALD","NEB","Direct")

MCBiasTab50 <- rbind(tapply(apply((w50JCs[1:maxcnt,] - qhats.pop5[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((q5bs[1:maxcnt,] - qhats.pop5[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((EB.norm.50s[1:maxcnt,] - qhats.pop5[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((exp(qhats.samp5[1:maxcnt,]) - qhats.pop5[1:maxcnt,]), 2, mean), nis, mean))
rownames(MCBiasTab50) <- c("LIGPD","ALD","NEB","Direct")

MCBiasTab75 <- rbind(tapply(apply((w75JCs[1:maxcnt,] - qhats.pop75[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((q75bs[1:maxcnt,] - qhats.pop75[1:maxcnt,]), 2, mean), nis, mean),
	tapply(apply((EB.norm.75s[1:maxcnt,] - qhats.pop75[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((exp(qhats.samp75[1:maxcnt,]) - qhats.pop75[1:maxcnt,]), 2, mean), nis, mean))
rownames(MCBiasTab75) <- c("LIGPD","ALD","NEB","Direct")

MCBiasTab90 <- rbind(tapply(apply((w90JCs[1:maxcnt,] - qhats.pop90[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((q90bs[1:maxcnt,] - qhats.pop90[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((EB.norm.90s[1:maxcnt,] - qhats.pop90[1:maxcnt,]) , 2, mean), nis, mean),
	tapply(apply((exp(qhats.samp90[1:maxcnt,]) - qhats.pop90[1:maxcnt,]), 2, mean), nis, mean))
rownames(MCBiasTab90) <- c("LIGPD","ALD","NEB","Direct")

##########  Output, including bootstrap bootstrap:

outmse1 <- cbind(tapply(apply( mhb25s, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply( mhb25s, 2, sd), nis[1:D[1]], mean)/sqrt(nrow(mhb25s)*20)*100,
tapply(apply( w25JCs  + sqrt(mhb25s)*1.96 > qhats.pop25  &  w25JCs  - sqrt(mhb25s)*1.96 < qhats.pop25 , 2, mean), nis[1:D[1]], mean),
tapply(apply((qhats.pop25  - w25JCs)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25  - EB.norm.25s)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25  - q25bs)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25  -  qhats.samp25)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25  - w25JCs) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25  - EB.norm.25s) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25  - q25bs) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop25  -  exp(qhats.samp25)) , 2, mean), nis[1:D[1]], mean)*100 

)



outmse2 <- cbind(tapply(apply( mhb50s, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply( mhb50s, 2, sd), nis[1:D[1]], mean)/sqrt(nrow(mhb25s)*20)*100,
tapply(apply( w50JCs  + sqrt(mhb50s)*1.96 > qhats.pop5  &  w50JCs  - sqrt(mhb50s)*1.96 < qhats.pop5 , 2, mean), nis[1:D[1]], mean),
tapply(apply((qhats.pop5  - w50JCs)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5  - EB.norm.50s)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5  - q5bs)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5  -  qhats.samp5)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5  - w50JCs) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5  - EB.norm.50s) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5  - q5bs) , 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop5  - exp(qhats.samp5)) , 2, mean), nis[1:D[1]], mean)*100 
)



 
outmse3 <- cbind(tapply(apply( mhb75s, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply( mhb75s, 2, sd), nis[1:D[1]], mean)/sqrt(nrow(mhb25s)*20)*100,
tapply(apply( w75JCs  + sqrt(mhb75s)*1.96 > qhats.pop75  &  w75JCs  - sqrt(mhb75s)*1.96 < qhats.pop75 , 2, mean), nis[1:D[1]], mean),
tapply(apply((qhats.pop75  - w75JCs)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75  - EB.norm.75s)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75  - q75bs)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75  -  qhats.samp75)^2, 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75  - w75JCs), 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75  - EB.norm.75s), 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75  - q75bs), 2, mean), nis[1:D[1]], mean)*100,
tapply(apply((qhats.pop75  -  exp(qhats.samp75)) , 2, mean), nis[1:D[1]], mean)*100 
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

