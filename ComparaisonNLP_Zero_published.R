set.seed(123456)

n = 10
sigma = 1




lBF_JR = function(x,n){ 
  .5*log(sigma) + .5*(x^2)*(n^2)/(1/sigma + n) - 0.5*log(1/sigma + n)
  log( (n*x)^2/(1/sigma + n)^2  + 1/(1/sigma + n) )
}

library(mombf)

lBF_RT = function(x,n){
  fitlm = lm(x~1)
  mombf(fitlm,coef = c(1),g = c(0.358),logbf = 1)
  
}


Stest = function(x,n){
  v = max((log(log(n))),1.96) #sqrt(log(n)/2)
  pnorm(-v/sqrt(n),m = n*x/(n+1/sigma),sd = 1/sqrt(n+1/sigma)) + 
    1-pnorm(v/sqrt(n),m = n*x/(n+1/sigma),sd = 1/sqrt(n+1/sigma))
}

lBF = function(x,n){ 
  -0.5*log(1+sigma*n) + (n*sigma*x^2)/(2*(sigma + 1/n))
}


test = function(n,theta){
  x = rnorm(1,m = theta,sd = 1/sqrt(n))
  return(c(lBF_JR(x,n)>0, Stest(x,n)>1/2,lBF(x,n) >0))
} 

nlist = c(100,250,500,10e3,5*10e3,10e4,5*10e4,10e5)
run = function(n) rowMeans(replicate(50000,test(n,sqrt(2*log(n)/n))))


run0 = function(n) rowMeans(replicate(50000,test(n,sqrt(log(n)/n))))


run1 = function(n) rowMeans(replicate(50000,test(n,0)))

m = t(sapply(nlist,run))



m2 = t(sapply(nlist,run0))

m3 = t(sapply(nlist,run1))

