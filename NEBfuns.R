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

cnp <- function(t, muvec, sdvec){
  pnorm(t, muvec, sdvec)
}





