# Title     : TODO
# Objective : TODO
# Created by: Emily
# Created on: 6/24/2017

betahat.glsfun <- function(sigma2vec, r.ols, k, X.mat, response){
    m <- length(r.ols)
    ci <- 1/sqrt(sigma2vec)/sum(1/sqrt(sigma2vec))
    sig2.tildeu <- sum(ci*(m/(m-k)*r.ols^2 - sigma2vec))
    Vhat1 <- sum(ci^2*((m/(m-k)*r.ols^2 - sigma2vec)-sig2.tildeu)^2)
    sig2u <- max(c(0.5*sqrt(Vhat1), sig2.tildeu))
    W.mat <- diag(1/(sig2u + sigma2vec))
    betahat.gls <- solve(t(X.mat)%*%W.mat%*%X.mat)%*%(t(X.mat)%*%W.mat%*%response)
    #sig2.tildeu
    list(betahat.gls, sig2u, sig2.tildeu)
}

