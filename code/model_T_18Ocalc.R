model {
  
  #Data model for d18O observations
  
  for(i in min(ind.up):max(ind.up)){
    d18OData[i] ~ dnorm(d18O_m[d18O.ind[i]], pre.d18O)
  }
  
  for(i in min(ind.mi):max(ind.mi)){
    d18OData[i] ~ dnorm(d18O_m[d18O.ind[i]], pre.d18O)
  }
  
  for(i in min(ind.lo):max(ind.lo)){
    d18OData[i] ~ dnorm(d18O_m[d18O.ind[i]], pre.d18O)
  }
  
  pre.d18O ~ dgamma(d18O.pre.sh, d18O.pre.ra)
  
  ###process model for benthic foram d18O###
  #with uncertainty assigned to each parameter
  #calculate d18O from temperature using equations from Hansen et al 2013
  
  for(i in min(tdo.lo):max(tdo.lo)){
    d18O_m[i] = (t_do[i] - a.m[1])/b.m[1] - c.m[1]
  }
  #age 3.600 to 34.025
  for(i in min(tdo.mi):max(tdo.mi)){
    d18O_m[i] = (t_do[i] - a.m[2])/b.m[2] - c.m[2]
  }
  #age 34.025 to 67.000
  for(i in min(tdo.up):max(tdo.up)){
    d18O_m[i] = (t_do[i] - a.m[3])/b.m[3] - c.m[3]
  }

  #modeling uncertainty in a, b and c
  a.m[1] ~ dnorm(a[1], 1/(sigma*a[1])^2)
  a.m[2] ~ dnorm(a[2], 1/(sigma*a[2])^2)
  a.m[3] ~ dnorm(a[3], 1/(sigma*a[3])^2)

  b.m[1] ~ dnorm(b[1], 1/(sigma*b[1])^2)
  b.m[2] ~ dnorm(b[2], 1/(sigma*b[2])^2)
  b.m[3] ~ dnorm(b[3], 1/(sigma*b[3])^2)

  c.m[1] ~ dnorm(c[1], pre.c)
  c.m[2] ~ dnorm(c[2], pre.c)
  c.m[3] ~ dnorm(c[3], pre.c)
  
  pre.c = 1/0.05^2 #c[3] = 0, so all cs are modeled to have the same uncertainty

  #according to Hansen et al. 2013:
  #age 0.000 to 3.660
  #t_do = 1 – 4.4 / 3 * (d18O_m – 3.25)
  
  #age 3.600 to 34.025
  #t_do = 5 – 8 / 3 * (d18O_m – 1.75)
  
  #age 34.025 to 67.000
  #t_do = 12 − 4 * (d18O_m + 0)
  #equation t_do = a + b * (d18O_m + c) according to Hansen et al. 2013
  #solve for d18O_m: d18O_m = (t_do - a)/b - c
  #assign a, b and c values
  
  a[1] = 1
  a[2] = 5
  a[3] = 12
  
  b[1] = -4.4/3
  b[2] = -8/3
  b[3] = -4
  
  c[1] = -3.25
  c[2] = -1.75
  c[3] = 0

  ###process model for deep ocean temperature###
  #surface ocean temp to deep ocean temp
  #ages 0.000 to 1.810
  for(i in min(ts.lo):max(ts.lo)){
    t_do[i] = (t_s[i] - y.m[1])/x.m[1]
  }
  #ages 1.810 to 5.330
  for(i in min(ts.mi):max(ts.mi)){
    t_do[i] = (t_s[i] - y.m[2])/x.m[2]
  }
  #ages 5.330 to 67.000
  for(i in min(ts.up):max(ts.up)){
    t_do[i] = (t_s[i] - y.m[3])/x.m[3]
  }

  #use the same sigma value here
  x.m[1] ~ dnorm(x[1], 1/(sigma*x[1])^2)
  x.m[2] ~ dnorm(x[2], 1/(sigma*x[2])^2)
  x.m[3] ~ dnorm(x[3], 1/(sigma*x[3])^2)
  
  y.m[1] ~ dnorm(y[1], 1/(sigma*y[1])^2)
  y.m[2] ~ dnorm(y[2], 1/(sigma*y[2])^2)
  y.m[3] ~ dnorm(y[3], 1/(sigma*y[3])^2)

  #according to Hansen et al. 2013:
  #ages 0.000 to 1.810
  #t_m = 2 * t_do + 12.25
  
  #ages 1.810 to 5.330
  #t_m = 2.5 * t_do + 12.15
  
  #ages 5.330 to 67.000
  #t_m = t_do + 14.15
  #equation: ts = x * td + y
  #solve for t_do: t_do = (t_m - y)/x
  #assign x and y values
  x[1] = 2
  x[2] = 2.5
  x[3] = 1
  
  y[1] = 12.25
  y[2] = 12.15
  y[3] = 14.15

  ###environmental model###
  #Convert to delta temperature GMST for comparison with the Westerhold stack
  GMST[1:al] = t_s[1:al] - ho.m

  ho.m ~ dnorm (ho.mean, 1/(sigma*ho.mean)^2)

  ho.mean = 14.15

  sigma = 0.01 #modeling 1 sd of parameters to be 1% of all parameters
  
  #temperature time series
  for(i in 2:al){
    t_s[i] = t_s[i-1] + t_s.eps[i]
    
    t_s.eps[i] ~ dnorm(t_s.eps[i-1] * t_s.eps.ac, t_s.pre) 
  }
  
  t_s.eps[1] ~ dnorm(0, t_s.pre)
  t_s[1] ~ dunif(15, 30)
  
  #Priors on model parameters  
  t_s.eps.ac ~ dunif(0.01, 0.99)
  
  t_s.pre ~ dgamma(t_s.pre.shp, t_s.pre.rate)
  t_s.pre.shp = 2
  t_s.pre.rate = 2
  
}