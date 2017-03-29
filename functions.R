g1 = function(x) - ((x<=1/2)*4*(x-0.5)**3.0 + 0.1*(x-0.5) - 0.25*exp(-250*(x-0.25)**2) )

g2= function(x)  x/10.0

g3= function(x) -( -0.1*exp( -50.0*(x-0.5)**2.0) )

g4= function(x)-(0.1*cos(6*pi*x)   )

g5= function(x) -( x/5 - g3(x) )

g6= function(x) -( x/5 - g4(x)  )

g7= function(x) -( x + 1 -0.25*exp(-50*(x - 0.5)**2)  )

g8= function(x) -x**2.0 / 2.0

g9= function(x) 0*x

g10= function(x) -x - 1

g11= function(x) -x -1 + (9.0/20)*exp(-50*(x-0.5)**2.0)


f1 = function(x){ 
  (-15*(x-0.5)**3)*(x <= 0.5)*1 - 0.3*(x - 0.5) + exp(-250*((x-0.25)**2))
}
f2= function(x){ 
  0.15*x 
}
f3= function(x){ 
  0.2*exp(-50*((x-0.5)**2))
}
f4= function(x){ 
  -0.1*cos(6*pi*x)*5 
}
f5= function(x){ 
  -0.2*x + f3(x) 
}
f6= function(x){ 
  -0.2*x + f4(x)
}
f7= function(x){ 
  -(1+x) + 0.45*exp(-50*((x-0.5)**2)) 
}
f02= function(x){ 
  -0.5*(x**2) 
}
f03= function(x){ 
  0.0*x
}

f04 = function(x){
  10*(x<0.5) 
}

f05 = function(x){
  -.75*x**2 + x**2*exp(-3*abs(x-0.5)^.5)
}

