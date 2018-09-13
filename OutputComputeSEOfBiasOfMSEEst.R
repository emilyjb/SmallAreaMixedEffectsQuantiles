
tstatrat <- function(num, den){
	denhat <- apply(den, 2, mean)
	rhat <- apply(num, 2, mean)/denhat
	devrat <- t(1/denhat*(t(num) -  rhat*t(den)))
	sehatrat <- apply(devrat, 2, sd)/sqrt(nrow(num))
	(rhat-1)/sehatrat
}


###  Approximate test statistic for bias of MSE estimator:

bhatmse25 <- apply( mhb25s - (qhats.pop25  - w25JCs)^2, 2, mean)
sebhatmse25 <- apply( mhb25s - (qhats.pop25  - w25JCs)^2, 2, sd)/sqrt(nrow(mhb25s))
sebhatmse25av <- sqrt(tapply(sebhatmse25^2, nis[1:D[1]], sum))/D[1]
bhatmse25av <- tapply(bhatmse25, nis[1:D[1]], mean)
tstatmse25 <- bhatmse25/sebhatmse25
tstatmse25av <- bhatmse25av/sebhatmse25av
tstatratmse25 <- tstatrat(mhb25s, (qhats.pop25  - w25JCs)^2)
mean(abs(tstatratmse25)>2)
boxplot(tstatmse25~nis[1:D[1]])
boxplot(tstatratmse25~nis[1:D[1]])

bhatmse50 <- apply( mhb50s - (qhats.pop5  - w50JCs)^2, 2, mean)
sebhatmse50 <- apply( mhb50s - (qhats.pop5  - w50JCs)^2, 2, sd)/sqrt(nrow(mhb50s))
sebhatmse50av <- sqrt(tapply(sebhatmse50^2, nis[1:D[1]], sum))/D[1]
bhatmse50av <- tapply(bhatmse50, nis[1:D[1]], mean)
tstatmse50 <- bhatmse50/sebhatmse50
tstatmse50av <- bhatmse50av/sebhatmse50av
tstatratmse50 <- tstatrat(mhb50s, (qhats.pop5  - w50JCs)^2)
mean(abs(tstatratmse50)>2)
boxplot(tstatmse50~nis[1:D[1]])
boxplot(tstatratmse50~nis[1:D[1]])



bhatmse75 <- apply( mhb75s - (qhats.pop75  - w75JCs)^2, 2, mean)
sebhatmse75 <- apply( mhb75s - (qhats.pop75  - w75JCs)^2, 2, sd)/sqrt(nrow(mhb75s))
sebhatmse75av <- sqrt(tapply(sebhatmse75^2, nis[1:D[1]], sum))/D[1]
bhatmse75av <- tapply(bhatmse75, nis[1:D[1]], mean)
tstatmse75 <- bhatmse75/sebhatmse75
tstatmse75av <- bhatmse75av/sebhatmse75av
tstatratmse75 <- tstatrat(mhb75s, (qhats.pop75  - w75JCs)^2)
mean(abs(tstatratmse75)>2)
boxplot(tstatmse75~nis[1:D[1]])
boxplot(tstatratmse75~nis[1:D[1]])


