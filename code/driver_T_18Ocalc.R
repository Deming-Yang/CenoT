#####
#Preliminaries
#####

#Load libraries
library(rjags)
library(R2jags)
library(openxlsx)
library(EnvStats)

#Read proxy data
df = "data/Westerhold.xlsx"
d = read.xlsx(df, sheet = "data")

#Data subset 
d = d[,c(1:3)]
names(d) = c("age","d18O", "temp")

#Set up ages vector
ages.bin = 0.5
ages = seq(70, 0, by = 0 - ages.bin) - ages.bin / 2
ages.len = length(ages)

#Age index
d$ai = ceiling((70 - d$age) / ages.bin)

#calculating 1 sigma of original data of each time bin
d18O.sd.tb <- rep(0,ages.len) #initiating vector
for(i in 1:ages.len){
  d18O.sd.tb[i] <- sd(d$d18O[which(d$ai == i)])
}
#this is used in the precition term of the data evaluation
hist(1/d18O.sd.tb^2)

egamma.d18O.sd.tb <- egamma(1/d18O.sd.tb^2)
d18O.pre.sh <- egamma.d18O.sd.tb$parameters[1]
d18O.pre.ra <- egamma.d18O.sd.tb$parameters[2]

#age index for deep ocean temp time series evaluation
#this is to account for the different equations used in the model calculation
tdo.lo <- which(ages < 3.660)
tdo.up <- which(ages > 34.025)
tdo.mi <- c((max(tdo.up)+1):(min(tdo.lo)-1))

#age index for surface ocean temp time series evaluation
#this is to account for the different equations used in the model calculation
ts.lo <- which(ages < 1.810)
ts.up <- which(ages > 5.330)
ts.mi <- c((max(ts.up)+1):(min(ts.lo)-1))

#age index for d18O evaluation age 0.000 to 3.660
ind.lo <- which(d$age < 3.660)

ind.up <- which(d$age > 34.025)

ind.mi <- c((max(ind.lo)+1):(min(ind.up)-1))


##Data to pass to BUGS model
dat = list("d18OData" = d$d18O, "d18O.ind" = d$ai, "al" = ages.len, 
           "ind.lo" = ind.lo, "ind.mi" = ind.mi, "ind.up" = ind.up,
           "tdo.lo" = tdo.lo, "tdo.mi" = tdo.mi, "tdo.up" = tdo.up,
           "ts.lo" = ts.lo, "ts.mi" = ts.mi, "ts.up" = ts.up,
           "d18O.pre.sh" = d18O.pre.sh, "d18O.pre.ra" = d18O.pre.ra)

##Parameters to save
parameters = c("t_s", "t_s.pre", "t_s.eps.ac","d18O_m","a.m", "b.m", "c.m", "x.m", "y.m", "GMST","ho.m")

##Run it
n.iter = 12000
n.burnin = 2000
n.thin = trunc((n.iter - n.burnin) / 2500)
pt = proc.time()
p_18Ocal = do.call(jags.parallel, list(model.file = "code/model_T_18Ocalc.R", parameters.to.save = parameters, 
                                      data = dat, inits = NULL, n.chains = 4, n.iter = n.iter, 
                                      n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt #~1 hour

save(p_18Ocal, file = "out/postd18OcalTemp.rda")

######use a smaller time bin, 0.2 Ma########
#Read proxy data
df = "data/Westerhold.xlsx"
d = read.xlsx(df, sheet = "data")

#Data subset 
d = d[,c(1:3)]
names(d) = c("age","d18O", "temp")

#Set up ages vector
ages.bin = 0.2
ages = seq(70, 0, by = 0 - ages.bin) - ages.bin / 2
ages.len = length(ages)

#Age index
d$ai = ceiling((70 - d$age) / ages.bin)

#calculating 1 sigma of original data of each time bin
d18O.sd.tb <- rep(0,ages.len) #initiating vector
for(i in 1:ages.len){
  d18O.sd.tb[i] <- sd(d$d18O[which(d$ai == i)])
}
#this is used in the precition term of the data evaluation
hist(1/d18O.sd.tb^2)

egamma.d18O.sd.tb <- egamma(1/d18O.sd.tb^2)
d18O.pre.sh <- egamma.d18O.sd.tb$parameters[1]
d18O.pre.ra <- egamma.d18O.sd.tb$parameters[2]

#age index for deep ocean temp time series evaluation
#this is to account for the different equations used in the model calculation
tdo.lo <- which(ages < 3.660)
tdo.up <- which(ages > 34.025)
tdo.mi <- c((max(tdo.up)+1):(min(tdo.lo)-1))

#age index for surface ocean temp time series evaluation
#this is to account for the different equations used in the model calculation
ts.lo <- which(ages < 1.810)
ts.up <- which(ages > 5.330)
ts.mi <- c((max(ts.up)+1):(min(ts.lo)-1))

#age index for d18O evaluation age 0.000 to 3.660
ind.lo <- which(d$age < 3.660)

ind.up <- which(d$age > 34.025)

ind.mi <- c((max(ind.lo)+1):(min(ind.up)-1))


##Data to pass to BUGS model
dat = list("d18OData" = d$d18O, "d18O.ind" = d$ai, "al" = ages.len, 
           "ind.lo" = ind.lo, "ind.mi" = ind.mi, "ind.up" = ind.up,
           "tdo.lo" = tdo.lo, "tdo.mi" = tdo.mi, "tdo.up" = tdo.up,
           "ts.lo" = ts.lo, "ts.mi" = ts.mi, "ts.up" = ts.up,
           "d18O.pre.sh" = d18O.pre.sh, "d18O.pre.ra" = d18O.pre.ra)

##Parameters to save
parameters = c("t_s", "t_s.pre", "t_s.eps.ac","d18O_m","a.m", "b.m", "c.m", "x.m", "y.m", "GMST","ho.m")

##Run it
n.iter = 12000
n.burnin = 2000
n.thin = trunc((n.iter - n.burnin) / 2500)
pt = proc.time()
p_18Ocal0.2 = do.call(jags.parallel, list(model.file = "code/model_T_18Ocalc.R", parameters.to.save = parameters, 
                                       data = dat, inits = NULL, n.chains = 4, n.iter = n.iter, 
                                       n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt #~2.5 hours

save(p_18Ocal0.2, file = "out/postd18OcalTemp0_2.rda")
