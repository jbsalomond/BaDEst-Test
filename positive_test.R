library(dplyr)
library(MASS)
library(invgamma)
library(pbapply)
library(parallel)
cl = cl <- makeCluster(10)
clusterExport(cl, ls())
clusterEvalQ(cl, library(MASS))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(invgamma))
clusterEvalQ(cl, library(pbapply))
# Function H 
H = function(omega){
 (- min(omega))
}

# sampling k 

## Computing the log prob 

Cxk = function(k,x,y,prior,alpha,beta,mu,lambda){
  n = length(x)
  TDC = mutate(data.frame(x), I = cut(x,breaks = seq(0,k)/k,right = F),Z = 1,Y = y)
  Yi2 = aggregate(Y ~ I, function(z){sum((z - mean(z))^2)},data = TDC)$Y
  #Yim = aggregate(Y ~ I, mean,data = TDC)$Y
  ni = aggregate(Z ~ I, sum,data = TDC)$Z
  #Y3 = ni*mu*(Yim - m)^2/(ni + mu)
  #beta = 2*mean(Yi)
  btilde = beta + 0.5*sum(Yi2) #+ 0.5*sum(Y3)
  if(prior == "pois"){ logprior = dpois(k-2,lambda,log = T) }
  if(prior == "geom"){ logprior = dgeom(k-1,lambda,log = T) }
  lpik = -log(btilde)*(alpha + n/2) - 0.5*sum(log(ni + mu)) + 0.5*k*log(mu/(2*pi)) +
    logprior
  return(lpik)
}

# sample k 
sampleK00 = function(x,y,prior,alpha,beta,mu,lambda,Kmax,Nk){
  kliste = seq(2,Kmax)
  logpik = sapply(X = kliste,FUN = Cxk,
                  x = x,y = y,prior = prior,alpha = alpha,beta = beta,mu = mu,lambda = lambda)
  pik = exp(logpik)/sum(exp(logpik))
  sample(2:Kmax,prob = pik,replace = T,size = Nk)
}


sampleK = function(x,y,prior,alpha,beta,mu,lambda,Kmax,Nk){
  k0 = ceiling(sqrt(length(x)))+2
  kliste = seq(2,k0)
  logpik = sapply(X = kliste,FUN = Cxk,
                  x = x,y = y,prior = prior,alpha = alpha,beta = beta,mu = mu,lambda = lambda)
  test = TRUE
  i = k0
  while(test){
    logpik = c(logpik,Cxk(k = i+1,x=x,y=y,prior = prior,alpha = alpha,beta = beta,mu=mu,lambda = lambda))
    #print(length(logpik)-k0)
    #print((logpik[i] - max(logpik)))
    test = (logpik[i] - max(logpik)) > -300
    #print(test)
    test = test&(i<Kmax)
    i = i+1
  }
  center = logpik - max(logpik)
  #print(center)
  pik = exp(center)/sum(exp(center))
  sample(2:i,prob = pik,replace = T,size = Nk)
}





test = function(x,y,prior,alpha,beta,mu,lambda,Kmax,Nk,M0,verbose=F,autoM0 = F){
  # Sample k first
  n = length(x)
  k = sampleK(x,y,prior,alpha,beta,mu,lambda,Kmax,Nk)
  df = data.frame(k = k, ni = 1)
  ktable = aggregate(ni~k,data = df,FUN = sum)
  #print(ktable)
  K = dim(ktable)[1]
  out = rep(0,K)
  if(verbose) print(ktable)
  for(i in 1:K){
    j = ktable[i,1]
    nj = ktable[i,2]
    outk = rep(0,nj)
    TDC = mutate(data.frame(x), I = cut(x,breaks = seq(0,j)/j,right = F),Z = 1,Y = y)
    Yi = aggregate(Y~I, data = TDC, mean)$Y
    ni = aggregate(Z~I, data = TDC, sum)$Z
    Yi2 = aggregate(Y ~ I, function(z){sum((z - mean(z))^2)},data = TDC)$Y
    #Yim = aggregate(Y ~ I, mean,data = TDC)$Y
    #Y3 = ni*mu*(Yim - m)^2/(ni + mu)
    btilde = beta + 0.5*sum(Yi2) #+ 0.5*sum(Y3)
    sigma = rinvgamma(nj,shape = alpha + n/2, rate = btilde)
    #print(median(sqrt(sigma)))
    if(verbose) print(median(sigma))
    
    for(l in 1:nj){
      if(autoM0) tau = 2*M0*sqrt(log(n/j))*(sqrt(j*median(sigma)/(n+j*mu)))
      else tau = (M0*sqrt(j*median(sigma)*log(n)/(n)))
      postmean = Yi
      postvar = sigma[l]/(ni + mu)
      omega = rnorm(n = j, m = postmean,sd = sqrt(postvar))
      #print(j)
      #print(omega)
      #print(H(omega))
      #print(tau)
      outk[l] = H(omega)>tau
    }
    out[i] = sum(outk)
  }
  return(sum(out)/Nk>=0.5)
}


run = function(n,sd = 0.1,verbose = T,autoM0 = F,prior = "geom",lambda = 0.4,M0=1,mu=1,Kmax=100){
sd = 0.1

ro = 6*(log(n)/n)^(1/3)
f = function(x){
  ro*(abs(x-0.5)-0.1)*(abs(x-0.5)<0.1)
}

X = seq(0,1,length = n+1)[-(n+1)]
y = f(X) + rnorm(n,sd = sd)
return(test(X,y,prior,alpha = 100, beta = .5,mu = mu,lambda = lambda,Kmax= Kmax,Nk = 2500,M0 = M0,verbose = verbose,autoM0 = autoM0) ) 
}

positive = function(spsize){
mean(replicate(run(n = spsize,
                   verbose = F,prior = "geom",lambda = .4,M0 = 1,mu = 1,autoM0 = T),n = 250))
}

positive(100)

nlist = c(100,250,500,1000,2500)


clusterExport(cl, ls())
start.time = proc.time()
resutl <- parLapply(cl,nlist,positive)
print(proc.time() - start.time)

t(unlist(resutl))