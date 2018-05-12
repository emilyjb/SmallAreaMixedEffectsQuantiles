bootgensimple2 <- function(sig2bhat.update, GN, XBhat.update1, CNis, nis, D, areafac.pop, rholxilrhouxiu){

  ahats.b <- rnorm(D, mean = 0, sd = sqrt(sig2bhat.update ))
  all.y.b <- XBhat.update1 + as.vector(GN%*%ahats.b)
  lyN.b <- areapredraninterp(all.y.b, rholxilrhouxiu, tauvec)	
  #ran.taus <-  sample(1:length(q.vec), size = nrow(all.y.b), replace = TRUE)
  #mu.boots <- sapply(1:nrow(all.y.b), function(i){ all.y.b[i,ran.taus[i]]})
  
#  lyN.b <- apply(all.y.b, 1, sample, size = 1, replace = FALSE)

  #t.min.boot <- min(all.y.b)
  #t.max.boot <- max(all.y.b)
  #t.seq.boot <- seq(t.min.boot, t.max.boot, length = length(q.vec))
  ####  Get the area quantile for the new all.y.b

  #q25b <- sapply(1:D, get.qpop.boot, 0.25, q.vec, all.y.b, t.min.boot, t.max.boot, t.seq.boot)
  #q50b <- sapply(1:D, get.qpop.boot, 0.50, q.vec, all.y.b, t.min.boot, t.max.boot, t.seq.boot)
  #q75b <- sapply(1:D, get.qpop.boot, 0.75, q.vec, all.y.b)

  #sapply(1:D, function(x){sample(1:Nx, size = )})

  smc.b <- sample(1:CNis[1], size=nis[1], replace=FALSE)
  i <- 1
  repeat{
    i <- i + 1
    smc.b <- c(smc.b, sample((CNis[(i-1)] +1):CNis[i], size=nis[i], replace=FALSE))
    if(i == D){break}
  }


  #lyN.b <- apply(all.y.b, 1, sample, size = 1)

  q25boot  <-  tapply(lyN.b, areafac.pop, quantile, prob = 0.25)
  q50boot <-   tapply(lyN.b, areafac.pop, quantile, prob = 0.5)
  q75boot <-  tapply(lyN.b, areafac.pop, quantile, prob = 0.75)

  #q75b <- c(q25b, q5b, q75b)

  lys.b <- lyN.b[smc.b]

  list(lys.b, smc.b, q25boot = q25boot, q50boot = q50boot, q75boot = q75boot)

}

bootgensimple2trans <- function(sig2bhat.update, GN, XBhat.update1, CNis, nis, D, areafac.pop, rholxilrhouxiu, lam){
  
  ahats.b <- rnorm(D, mean = 0, sd = sqrt(sig2bhat.update ))
  all.y.b <- XBhat.update1 + as.vector(GN%*%ahats.b)
  lyN.b <- areapredraninterp(all.y.b, rholxilrhouxiu, tauvec)	
  #ran.taus <-  sample(1:length(q.vec), size = nrow(all.y.b), replace = TRUE)
  #mu.boots <- sapply(1:nrow(all.y.b), function(i){ all.y.b[i,ran.taus[i]]})
  
  #  lyN.b <- apply(all.y.b, 1, sample, size = 1, replace = FALSE)
  
  #t.min.boot <- min(all.y.b)
  #t.max.boot <- max(all.y.b)
  #t.seq.boot <- seq(t.min.boot, t.max.boot, length = length(q.vec))
  ####  Get the area quantile for the new all.y.b
  
  #q25b <- sapply(1:D, get.qpop.boot, 0.25, q.vec, all.y.b, t.min.boot, t.max.boot, t.seq.boot)
  #q50b <- sapply(1:D, get.qpop.boot, 0.50, q.vec, all.y.b, t.min.boot, t.max.boot, t.seq.boot)
  #q75b <- sapply(1:D, get.qpop.boot, 0.75, q.vec, all.y.b)
  
  #sapply(1:D, function(x){sample(1:Nx, size = )})
  
  # smc.b <- sample(1:CNis[1], size=nis[1], replace=FALSE)
  # i <- 1
  # repeat{
  #  i <- i + 1
  #  smc.b <- c(smc.b, sample((CNis[(i-1)] +1):CNis[i], size=nis[i], replace=FALSE))
  #  if(i == D){break}
  #}
  
  
  #lyN.b <- apply(all.y.b, 1, sample, size = 1)
  #htransGJ1()
  q25boot  <-  tapply(hinvtransG1(lyN.b, lam), areafac.pop, quantile, prob = 0.25)
  q50boot <-   tapply(hinvtransG1(lyN.b, lam), areafac.pop, quantile, prob = 0.5)
  q75boot <-  tapply(hinvtransG1(lyN.b, lam), areafac.pop, quantile, prob = 0.75)
  
  #q75b <- c(q25b, q5b, q75b)
  
  lys.b <- lyN.b[smc]
  
  list(lys.b, smc, q25boot = q25boot, q50boot = q50boot, q75boot = q75boot)
  
}



bootgensimple2transstep <- function(sig2bhat.update, GN, XBhat.update1, CNis, nis, D, areafac.pop, rholxilrhouxiu, lam){
  
  ahats.b <- rnorm(D, mean = 0, sd = sqrt(sig2bhat.update ))
  all.y.b <- XBhat.update1 + as.vector(GN%*%ahats.b)
  #lyN.b <- areapredraninterp(all.y.b, rholxilrhouxiu, tauvec)	
  #ran.taus <-  sample(1:length(q.vec), size = nrow(all.y.b), replace = TRUE)
  #mu.boots <- sapply(1:nrow(all.y.b), function(i){ all.y.b[i,ran.taus[i]]})
  
  lyN.b <- apply(all.y.b, 1, sample, size = 1, replace = FALSE)
  
  #t.min.boot <- min(all.y.b)
  #t.max.boot <- max(all.y.b)
  #t.seq.boot <- seq(t.min.boot, t.max.boot, length = length(q.vec))
  ####  Get the area quantile for the new all.y.b
  
  #q25b <- sapply(1:D, get.qpop.boot, 0.25, q.vec, all.y.b, t.min.boot, t.max.boot, t.seq.boot)
  #q50b <- sapply(1:D, get.qpop.boot, 0.50, q.vec, all.y.b, t.min.boot, t.max.boot, t.seq.boot)
  #q75b <- sapply(1:D, get.qpop.boot, 0.75, q.vec, all.y.b)
  
  #sapply(1:D, function(x){sample(1:Nx, size = )})
  
 # smc.b <- sample(1:CNis[1], size=nis[1], replace=FALSE)
 # i <- 1
 # repeat{
 #   i <- i + 1
 #   smc.b <- c(smc.b, sample((CNis[(i-1)] +1):CNis[i], size=nis[i], replace=FALSE))
 # if(i == D){break}
 # }
  
  
  #lyN.b <- apply(all.y.b, 1, sample, size = 1)
  #htransGJ1()
  q25boot  <-  tapply(hinvtransG1(lyN.b, lam), areafac.pop, quantile, prob = 0.25)
  q50boot <-   tapply(hinvtransG1(lyN.b, lam), areafac.pop, quantile, prob = 0.5)
  q75boot <-  tapply(hinvtransG1(lyN.b, lam), areafac.pop, quantile, prob = 0.75)
  
  #q75b <- c(q25b, q5b, q75b)
  
  lys.b <- lyN.b[smc]
  
  list(lys.b, smc, q25boot = q25boot, q50boot = q50boot, q75boot = q75boot)
  
}



bootgensimple3 <- function(sig2bhat.update, GN, XBhat.update1, CNis, nis, D, areafac.pop, rholxilrhouxiu, smc, b.dist){
if(b.dist == "Normal"){
  ahats.b <- rnorm(D, mean = 0, sd = sqrt(sig2bhat.update ))
}else{
  Ul <- runif(D, -0.5, 0.5)
  b <- sqrt(sig2bhat.update/2)
  ahats.b <- -b*sign(Ul)*log(1 - 2*abs(Ul))
}
  all.y.b <- XBhat.update1 + as.vector(GN%*%ahats.b)
  lyN.b <- areapredraninterp(all.y.b, rholxilrhouxiu, tauvec)	
  #ran.taus <-  sample(1:length(q.vec), size = nrow(all.y.b), replace = TRUE)
  #mu.boots <- sapply(1:nrow(all.y.b), function(i){ all.y.b[i,ran.taus[i]]})
  
#  lyN.b <- apply(all.y.b, 1, sample, size = 1, replace = FALSE)

  #t.min.boot <- min(all.y.b)
  #t.max.boot <- max(all.y.b)
  #t.seq.boot <- seq(t.min.boot, t.max.boot, length = length(q.vec))
  ####  Get the area quantile for the new all.y.b

  #q25b <- sapply(1:D, get.qpop.boot, 0.25, q.vec, all.y.b, t.min.boot, t.max.boot, t.seq.boot)
  #q50b <- sapply(1:D, get.qpop.boot, 0.50, q.vec, all.y.b, t.min.boot, t.max.boot, t.seq.boot)
  #q75b <- sapply(1:D, get.qpop.boot, 0.75, q.vec, all.y.b)

  #sapply(1:D, function(x){sample(1:Nx, size = )})

  #smc.b <- sample(1:CNis[1], size=nis[1], replace=FALSE)
  #i <- 1
  #repeat{
 #   i <- i + 1
 #   smc.b <- c(smc.b, sample((CNis[(i-1)] +1):CNis[i], size=nis[i], replace=FALSE))
 #   if(i == D){break}
 # }


  #lyN.b <- apply(all.y.b, 1, sample, size = 1)

  q25boot  <-  tapply(lyN.b, areafac.pop, quantile, prob = 0.25)
  q50boot <-   tapply(lyN.b, areafac.pop, quantile, prob = 0.5)
  q75boot <-  tapply(lyN.b, areafac.pop, quantile, prob = 0.75)

  #q75b <- c(q25b, q5b, q75b)

  lys.b <- lyN.b[smc]

  list(lys.b, smc, q25boot = q25boot, q50boot = q50boot, q75boot = q75boot)

}


