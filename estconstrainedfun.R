

estconstrained <- function(dat.temp, fml, tauvec){
	rqCur <- rq(fml, data = dat.temp, tau = 0.5, method = "fn")
	coefCur <- coef(rqCur)
	fitCur <- fitted(rqCur)
	iter.taus <- which(tauvec == 0.5)
	coefalls <- coefCur
	repeat{
		iter.taus <- iter.taus + 1
		rqUpdate <- rq(fml, data = dat.temp, tau = tauvec[iter.taus ], method = "fnc", R = lXN, r = lXN%*%coefCur)
		coefCur <- coef(rqUpdate)
		coefalls <- cbind(coefalls, coefCur)
		rqCur <- rqUpdate
		print(iter.taus)
		if(iter.taus == length(tauvec)){break}
	}
	rqCur <- rq(fml, data = dat.temp, tau = 0.5, method = "fn")
	coefCur <- coef(rqCur)
	fitCur <- fitted(rqCur)
	iter.taus <- which(tauvec == 0.5)
	repeat{
		iter.taus <- iter.taus - 1
		rqUpdate <- rq(fml, data = dat.temp, tau = tauvec[iter.taus ], method = "fnc", R = -lXN, r = -lXN%*%coefCur)
		coefCur <- coef(rqUpdate)
		coefalls <- cbind( coefCur, coefalls)
		rqCur <- rqUpdate
		print(iter.taus)
		if(iter.taus == 1){break}
	}
	coefalls
}



