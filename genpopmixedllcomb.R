

####  genpopmixedllf

genpopmixedll.comb.ebdist <- function(D, CNis, sig2le, sig2lu, mull0, beta1, Nis, lxN, GN, dfpick, type, bdist, tau){
   N <- CNis[D]
  if(type == "T"){
    leN <- (rt(N, df = dfpick) )/sqrt(dfpick/(dfpick - 2))
  }
  if(type == "Chi"){
    leN <- (rchisq(N, df = dfpick)-dfpick)/sqrt(2*dfpick)*(1+0.1*lxN)
  }
if(type == "SN"){
      alpha <- -5
    delta <- alpha/sqrt(1+alpha^2)
     w <- sqrt(1/(1-2*delta^2/pi))
      xi <- -w*delta*sqrt(2/pi)
	leN <- rsn(N,  xi, w, alpha)
 }
if(type == "Mix"){
      indicatbigvar <- rbinom(1, N, 0.2)
	leN <- rnorm(N, 0, 1)*(1-indicatbigvar) + rnorm(N, 0, 2)*indicatbigvar
}
   if(bdist == "Laplace"){
     Ul <- runif(D, -0.5, 0.5)
     b <- sqrt(sig2lu/2)
     uis <- -b*sign(Ul)*log(1 - 2*abs(Ul))
   }
   if(bdist == "Normal"){
    uis <- rnorm(D)*sqrt(sig2lu) 
   }
  luN <- GN%*%(uis)
  lyN <- mull0 + beta1*lxN + luN + leN*sqrt(sig2le)
  YN <- exp(lyN)
  YbarNis <- t(GN)%*%YN/Nis
  qij <-   ( mull0 + beta1*lxN + luN) + qt(tau, df=5)
  s <- sample(1:CNis[1], size=nis[1], replace=FALSE)
  i <- 1
  repeat{
    i <- i + 1
    s <- c(s, sample((CNis[(i-1)] +1):CNis[i], size=nis[i], replace=FALSE))
    if(i == D){break}
  }
  list(lyN, YN, YbarNis, s, qij, uis)
}

genpopmixedll <- function(D, CNis, sig2le, sig2lu, mull0, beta1, Nis, lxN, GN, tau = 0.5){
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
	  qij <- qnorm(tau, mean = ( mull0 + beta1*lxN + luN), sd = sqrt(sig2le))
        list(lyN, YN, YbarNis, s, qij)
}


genpopmixedll.t <- function(D, CNis, sig2le, sig2lu, mull0, beta1, Nis, lxN, GN, tau = 0.5){
        N <- CNis[D]
        leN <- rt(N, df = 5)/sqrt(5/3) 
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
	  qij <-   ( mull0 + beta1*lxN + luN) + qt(tau, df=1000)
        list(lyN, YN, YbarNis, s, qij)
}





genpopmixedll.comb <- function(D, CNis, sig2le, sig2lu, mull0, beta1, Nis, lxN, GN, tau = 0.5, dfpick, type){
	  
        N <- CNis[D]
	  if(type == "T"){
	        leN <- (rt(N, df = dfpick) )/sqrt(dfpick/(dfpick - 2))
	  }
	  if(type == "Chi"){
	        leN <- (rchisq(N, df = dfpick)-dfpick)/sqrt(2*dfpick)
	  }
	  uis <- rnorm(D)*sqrt(sig2lu) 
        luN <- GN%*%(uis)
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
        qij <-   ( mull0 + beta1*lxN + luN) + qt(tau, df=5)
        list(lyN, YN, YbarNis, s, qij, uis)
}


genpopmixedll.snc.comb <- function(D, CNis, sig2le, sig2lu, mull0, beta1, Nis, lxN, GN, tau = 0.5, delta,  dfpick, type){
        N <- CNis[D]
	  if(type == "T"){
	        leN <- (rt(N, df = dfpick))/sqrt(dfpick/(dfpick - 2))
	  }
	  if(type == "Chi"){
	        leN <- (rchisq(N, df = dfpick)-dfpick)/sqrt(2*dfpick)
	  }
	  W <- abs(rnorm(D))
	  omega2 <- sig2lu/(1-2/pi*delta^2)
	  alpha <- sqrt(delta^2/(1-delta^2))
	  xisn <- -sqrt(omega2)*delta*sqrt(2/pi)
	  uis <- rnorm(D, mean =  xisn + sqrt(omega2)*delta*W, sd = sqrt(omega2*(1-delta^2)))
#  Ul <- runif(D, -0.5, 0.5)
#	  b <- sqrt(sig2lu/2)
#	  uis <- -b*sign(Ul)*log(1 - 2*abs(Ul))
	  #uis <- rnorm(D)*sqrt(sig2lu) 
        luN <- GN%*%(uis)
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
  	qij <- (mull0 + beta1*lxN + luN) + qt(tau, df=5)
        list(lyN, YN, YbarNis, s, qij, uis)
}

#Ul <- runif(1000, -0.5, 0.5)
#b <- sqrt(0.36/2)
#-b*sign(Ul)*log(1- 2*abs(Ul))


genpopmixedll.lc.comb <- function(D, CNis, sig2le, sig2lu, mull0, beta1, Nis, lxN, GN, tau = 0.5, dfpick, type){
        N <- CNis[D]
#       leN <- (rchisq(N, df = 2) - 2)/sqrt(4)
	  if(type == "T"){
	        leN <- (rt(N, df = dfpick))/sqrt(dfpick/(dfpick - 2))
	  }
	  if(type == "Chi"){
	        leN <- (rchisq(N, df = dfpick)-dfpick)/sqrt(2*dfpick)
	  }
	  Ul <- runif(D, -0.5, 0.5)
	  b <- sqrt(sig2lu/2)
	  uis <- -b*sign(Ul)*log(1 - 2*abs(Ul))
	  #uis <- rnorm(D)*sqrt(sig2lu) 
        luN <- GN%*%(uis)
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
  qij <-   ( mull0 + beta1*lxN + luN) + qt(tau, df=5)
        list(lyN, YN, YbarNis, s, qij, uis)
}

#Ul <- runif(1000, -0.5, 0.5)
#b <- sqrt(0.36/2)
#-b*sign(Ul)*log(1- 2*abs(Ul))


genpopmixedll.stc.comb <- function(D, CNis, sig2le, sig2lu, mull0, beta1, Nis, lxN, GN, tau = 0.5, dfpick, type){
        N <- CNis[D]
#       leN <- (rchisq(N, df = 2) - 2)/sqrt(4)
	  if(type == "T"){
	        leN <- (rt(N, df = dfpick))/sqrt(dfpick/(dfpick - 2))
	  }
	  if(type == "Chi"){
	        leN <- (rchisq(N, df = dfpick)-dfpick)/sqrt(2*dfpick)
	  }
        uis <- rt.scaled(D, df = 3, mean = 0, sd = sqrt(sig2lu*(1)/3))
	 # Ul <- runif(D, -0.5, 0.5)
	 # b <- sqrt(sig2lu/2)
	 # uis <- -b*sign(Ul)*log(1 - 2*abs(Ul))
	  #uis <- rnorm(D)*sqrt(sig2lu) 
        luN <- GN%*%(uis)
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
  qij <-   ( mull0 + beta1*lxN + luN) + qt(tau, df=5)
        list(lyN, YN, YbarNis, s, qij, uis)
}



