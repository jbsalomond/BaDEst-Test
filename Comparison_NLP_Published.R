set.seed(123456)




Stest = function(x,n){
  sigma2 = 5
  v = max((log(log(n))),1.64)
  sp = sqrt(1/(n + 1/sigma2)) 
  mp = n*x*sp^2
  (1 - pnorm(v/sqrt(n),m = mp,sd = sp))>=1/2
  
}

BF = function(x,n){
  sigma2 =2
  sp = sqrt(1/(n + 1/sigma2)) 
  mp = n*x*sp^2
  BF = pnorm(-mp/sp,log = TRUE) - pnorm(mp/sp,log = TRUE)
  BF < 0 
  
}


JR = function(x,n){
  sigma2 = 0.358^2
  sp = sqrt(1/(n + 1/sigma2)) 
  mp = n*x*sp^2
  
  num = pnorm(0,m = mp,sd = sp,log = TRUE)
  denom = dnorm(mp/sp)*mp*sp + pnorm(mp/sp)*(mp^2+sp^2)
  (num - log(denom))<1
}


test = function(n,theta){
  x = rnorm(1,m = theta,sd = 1/sqrt(n))
  return(c(JR(x,n), Stest(x,n),BF(x,n)))
} 



run = function(n) rowMeans(replicate(50000,test(n,sqrt(2*log(n)/n))))


run0 = function(n) rowMeans(replicate(50000,test(n,sqrt(log(n)/n))))

run1 = function(n) rowMeans(replicate(50000,test(n,0)))


m = t(sapply(nlist,run))


m2 = t(sapply(nlist,run0))


m3 = t(sapply(nlist,run1))


