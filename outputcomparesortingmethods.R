
#########  Constrain
tapply(apply( (w10JCs - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JCs - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JCs - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JCs - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JCs - qhats.pop90)^2, 2, mean), nis, mean)

tapply(apply( (w10JC2s - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JC2s - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JC2s - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JC2s - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JC2s - qhats.pop90)^2, 2, mean), nis, mean)



#########  Sort
tapply(apply( (w10JSs - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JSs - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JSs - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JSs - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JSs - qhats.pop90)^2, 2, mean), nis, mean)

tapply(apply( (w10JS2s - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JS2s - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JS2s - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JS2s - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JS2s - qhats.pop90)^2, 2, mean), nis, mean)

tapply(apply( (w10JSs - qhats.pop10)^2, 2, mean), nis, mean)/tapply(apply( (w10JS2s - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JSs - qhats.pop25)^2, 2, mean), nis, mean)/tapply(apply( (w25JS2s - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JSs - qhats.pop5)^2, 2, mean), nis, mean)/tapply(apply( (w50JS2s - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JSs - qhats.pop75)^2, 2, mean), nis, mean)/tapply(apply( (w75JS2s - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JSs - qhats.pop90)^2, 2, mean), nis, mean)/tapply(apply( (w90JS2s - qhats.pop90)^2, 2, mean), nis, mean)


#########  Isoreg
tapply(apply( (w10JRs - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JRs - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JRs - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JRs - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JRs - qhats.pop90)^2, 2, mean), nis, mean)

tapply(apply( (w10JR2s - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JR2s - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JR2s - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JR2s - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JR2s - qhats.pop90)^2, 2, mean), nis, mean)

##### MSE ratios: constrain compare to JC2 in denominator:

tapply(apply( (w10JCs - qhats.pop10)^2, 2, mean), nis, mean)/tapply(apply( (w10JC2s - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JCs - qhats.pop25)^2, 2, mean), nis, mean)/tapply(apply( (w25JC2s - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JCs - qhats.pop5)^2, 2, mean), nis, mean)/tapply(apply( (w50JC2s - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JCs - qhats.pop75)^2, 2, mean), nis, mean)/tapply(apply( (w75JC2s - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JCs - qhats.pop90)^2, 2, mean), nis, mean)/tapply(apply( (w90JC2s - qhats.pop90)^2, 2, mean), nis, mean)

tapply(apply( (w10JS2s - qhats.pop10)^2, 2, mean), nis, mean)/tapply(apply( (w10JC2s - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JS2s - qhats.pop25)^2, 2, mean), nis, mean)/tapply(apply( (w25JC2s - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JS2s - qhats.pop5)^2, 2, mean), nis, mean)/tapply(apply( (w50JC2s - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JS2s - qhats.pop75)^2, 2, mean), nis, mean)/tapply(apply( (w75JC2s - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JS2s - qhats.pop90)^2, 2, mean), nis, mean)/tapply(apply( (w90JC2s - qhats.pop90)^2, 2, mean), nis, mean)


tapply(apply( (w10JR2s - qhats.pop10)^2, 2, mean), nis, mean)/tapply(apply( (w10JC2s - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JR2s - qhats.pop25)^2, 2, mean), nis, mean)/tapply(apply( (w25JC2s - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JR2s - qhats.pop5)^2, 2, mean), nis, mean)/tapply(apply( (w50JC2s - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JR2s - qhats.pop75)^2, 2, mean), nis, mean)/tapply(apply( (w75JC2s - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JR2s - qhats.pop90)^2, 2, mean), nis, mean)/tapply(apply( (w90JC2s - qhats.pop90)^2, 2, mean), nis, mean)



##### MSE ratios: check if isoreg or sort is better approximation to the constraining 1-step

tapply(apply( (w10JRs - qhats.pop10)^2, 2, mean), nis, mean)/tapply(apply( (w10JCs - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JRs - qhats.pop25)^2, 2, mean), nis, mean)/tapply(apply( (w25JCs - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JRs - qhats.pop5)^2, 2, mean), nis, mean)/tapply(apply( (w50JCs - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JRs - qhats.pop75)^2, 2, mean), nis, mean)/tapply(apply( (w75JCs - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JRs - qhats.pop90)^2, 2, mean), nis, mean)/tapply(apply( (w90JCs - qhats.pop90)^2, 2, mean), nis, mean)


tapply(apply( (w10JSs - qhats.pop10)^2, 2, mean), nis, mean)/tapply(apply( (w10JCs - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JSs - qhats.pop25)^2, 2, mean), nis, mean)/tapply(apply( (w25JCs - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JSs - qhats.pop5)^2, 2, mean), nis, mean)/tapply(apply( (w50JCs - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JSs - qhats.pop75)^2, 2, mean), nis, mean)/tapply(apply( (w75JCs - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JSs - qhats.pop90)^2, 2, mean), nis, mean)/tapply(apply( (w90JCs - qhats.pop90)^2, 2, mean), nis, mean)

tapply(apply( (w10JR2s - qhats.pop10)^2, 2, mean), nis, mean)/tapply(apply( (w10JC2s - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JR2s - qhats.pop25)^2, 2, mean), nis, mean)/tapply(apply( (w25JC2s - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JR2s - qhats.pop5)^2, 2, mean), nis, mean)/tapply(apply( (w50JC2s - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JR2s - qhats.pop75)^2, 2, mean), nis, mean)/tapply(apply( (w75JC2s - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JR2s - qhats.pop90)^2, 2, mean), nis, mean)/tapply(apply( (w90JC2s - qhats.pop90)^2, 2, mean), nis, mean)


tapply(apply( (w10JS2s - qhats.pop10)^2, 2, mean), nis, mean)/tapply(apply( (w10JC2s - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JS2s - qhats.pop25)^2, 2, mean), nis, mean)/tapply(apply( (w25JC2s - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JS2s - qhats.pop5)^2, 2, mean), nis, mean)/tapply(apply( (w50JC2s - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JS2s - qhats.pop75)^2, 2, mean), nis, mean)/tapply(apply( (w75JC2s - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JS2s - qhats.pop90)^2, 2, mean), nis, mean)/tapply(apply( (w90JC2s - qhats.pop90)^2, 2, mean), nis, mean)





##### MSE ratios: compare isoreg to sort 

tapply(apply( (w10JS2s - qhats.pop10)^2, 2, mean), nis, mean)/tapply(apply( (w10JR2s - qhats.pop10)^2, 2, mean), nis, mean)
tapply(apply( (w25JS2s - qhats.pop25)^2, 2, mean), nis, mean)/tapply(apply( (w25JR2s - qhats.pop25)^2, 2, mean), nis, mean)
tapply(apply( (w50JS2s - qhats.pop5)^2, 2, mean), nis, mean)/tapply(apply( (w50JR2s - qhats.pop5)^2, 2, mean), nis, mean)
tapply(apply( (w75JS2s - qhats.pop75)^2, 2, mean), nis, mean)/tapply(apply( (w75JR2s - qhats.pop75)^2, 2, mean), nis, mean)
tapply(apply( (w90JS2s - qhats.pop90)^2, 2, mean), nis, mean)/tapply(apply( (w90JR2s - qhats.pop90)^2, 2, mean), nis, mean)









