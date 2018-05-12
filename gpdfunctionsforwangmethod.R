



GInversePareto <- function(t, xi, rho){
  ( (1-t)^-xi -1)*rho/xi
}

GPareto <- function(y, xi, rho){
  1 - (1 + xi*y/rho)^(-(1/xi))
}

gpd.dens0 <- function(d, rho, xi, taul, tauu){
	if(xi != 0){
		gpdd <- (1 + xi*d/rho)^(-(1 + 1/xi))/rho
		if(xi > 0){
			gpdd <- gpdd*ifelse(d > -rho/xi, 1, 0)
		}
		if(xi < 0){
			gpdd <- gpdd*ifelse(d < -rho/xi, 1, 0)
		}
	
	}else{
		gpdd <- 1/rho*exp(-d/rho)
	}
	(gpdd/rho)
}

gpdxi <- function(xi, dvec, rho, taul, tauu){
	prod(sapply(dvec, gpd.dens,  rho, xi, taul, tauu))
}
	
gpd.dens <- function(d, rho, xi, taul, tauu){
	dgpd(d, loc = 0, scale= rho, shape = xi)
}

