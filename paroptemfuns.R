


opt.par.tau.em3Par <- function(tau, coef.start.all, bcond.ran, tauvec, lys, lXs, Gs){
        ind <- which(tauvec == tau)
        coef.start <- coef.start.all[,ind]
        optim(coef.start, condmean.ploss3Par, bcond.ran=bcond.ran, tau=tau, lys=lys, lXs=lXs, Gs=Gs)$par
}
 condmean.ploss3Par <- function(paravec, bcond.ran, tau, lys, lXs, Gs){ 
        mean(parApply(cl, bcond.ran, 2, ploss3, paravec, tau, lys, lXs, Gs))
}

