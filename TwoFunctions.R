### Generate lognormal data
genpopmixedll <- function(D, CNis, sig2le, sig2lu, mull0, beta1, Nis, lxN, GN){
	N <- CNis[D]
	leN <- rnorm(N)*sqrt(sig2le)
	luN <- GN%*%(rnorm(D)*sqrt(sig2lu))
	lyN <- mull0 + beta1*lxN + luN + leN
	YN <- exp(lyN)
	YbarNis <- t(GN)%*%YN/Nis
	s <- sample(1:CNis[1], size=nis[1], replace=FALSE)
	i <- 1
	repeat{
		i <- i + 1
		s <- c(s, sample((CNis[(i-1)] +1):CNis[i], size=nis[i], replace=FALSE))
		if(i == D){break}
	}
	list(lyN, YN, YbarNis, s)
}

### Compute REML estimates
getremlests <- function(s, D, YN, XN, GN){
        Ys <- YN[s]
        Xs <- XN[s,]
        Gs <- GN[s,]
        areagroup <- factor(Gs%*%(1:D))
        remlestsr <- lme(Ys~Xs[,2], random=~1|areagroup)
	  VC <- VarCorr(remlestsr)
        sig2uhat <- as.numeric(VC[1,1])
        sig2ehat <- as.numeric(VC[2,1])
        betahat <- matrix( (summary(remlestsr))$coef$fixed,2,1)
        vhatbetahat <- as.matrix(summary(remlestsr)$varFix)
        list(sig2uhat, sig2ehat, betahat, vhatbetahat)
}


#Repeat the elements of v1[i] v2[i] times
expandfun <- function(v1, v2){
v <- NULL
cnt <- 0
repeat{
cnt <- cnt +1
v <- c(v, rep(v1[cnt], times=v2[cnt]))
if(cnt==length(v1)){break}
}
v
}



vcovraneff2 <- function(x, sig2u, sig2e){
        ni <- sum(x)
        alphai <- sig2e + ni*sig2u
        Ivv <- ni^2/alphai^2
        Iee <- (ni-1)/sig2e^2 + 1/alphai^2
        Ive <- ni/alphai^2
        c(Ivv, Ive, Ive, Iee)
}




