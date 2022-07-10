#####
#Preliminaries
#####

##Environmental model for GMST time series, using the GMST data from the Westerhold stack
#Read proxy data
df = "data/Westerhold.xlsx"
d = read.xlsx(df, sheet = "data")

#Data subset 
d = d[,c(1:3)]
names(d) = c("age","d18O", "temp")
#make sure to load the "tsdens.R" function
ages.bin = 0.5
ages = seq(70, 0, by = 0 - ages.bin) - ages.bin / 2
ages.len = length(ages)

#compared to GMST stack by Westerhold
load("out/postd18OcalTemp.rda") #make sure to run "driver_T_18Ocalc.R" first
#compared to GMST stack by Westerhold
tp = p_18Ocal$BUGSoutput$sims.list$GMST
tp = tp[,-(ncol(tp))]
ages = ages[-length(ages)]#revert this vector

tpts = apply(tp, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
tpts = t(tpts)

plot(-10, 0, xlab="Age (Ma)", ylab = "GMST (relative to preindustrial)", 
     xlim=c(70,0), ylim = c(-5,25), axes = FALSE)

points(d$age, d$temp, cex=0.5, col = "dark grey")

tsdens(cbind(ages, tpts), "firebrick4")
axis(1)
axis(2)

#compared to corrected d18O benthic stack by Westerhold
#the d18O data is used to evaluate the proposed temp time series given the temp_18O equations
d18Op = p_18Ocal$BUGSoutput$sims.list$d18O_m
d18Op = d18Op[,-(ncol(d18Op))]

d18Opts = apply(d18Op, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
d18Opts = t(d18Opts)

plot(-10, 0, xlab="Age (Ma)", ylab = "d18O (‰, VSMOW)", 
     xlim=c(70,0), ylim = rev(range(d$d18O)), axes = FALSE)

points(d$age, d$d18O, cex=0.5, col = "grey")

tsdens(cbind(ages, d18Opts), "blue")
axis(1)
axis(2)

##########use a smaller time bin, 0.2 Ma#########
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

#compared to GMST stack by Westerhold
load("out/postd18OcalTemp0_2.rda") #make sure to run "driver_T_18Ocalc.R" first

#compared to GMST stack by Westerhold
tp0.2 = p_18Ocal0.2$BUGSoutput$sims.list$GMST
tp0.2 = tp0.2[,-(ncol(tp0.2))]
ages = ages[-length(ages)]#revert this vector

tpts0.2 = apply(tp0.2, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
tpts0.2 = t(tpts0.2)

plot(-10, 0, xlab="Age (Ma)", ylab = "GMST (relative to preindustrial)", 
     xlim=c(70,0), ylim = c(-5,25), axes = FALSE)

points(d$age, d$temp, cex=0.5, col = "dark grey")

tsdens(cbind(ages, tpts0.2), "firebrick4")
axis(1)
axis(2)

#compared to corrected d18O benthic stack by Westerhold
d18Op0.2 = p_18Ocal0.2$BUGSoutput$sims.list$d18O_m
d18Op0.2 = d18Op0.2[,-(ncol(d18Op0.2))]

d18Opts0.2 = apply(d18Op0.2, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
d18Opts0.2 = t(d18Opts0.2)

plot(-10, 0, xlab="Age (Ma)", ylab = "d18O (‰, VSMOW)", 
     xlim=c(70,0), ylim = rev(range(d$d18O)), axes = FALSE)

points(d$age, d$d18O, cex=0.5, col = "grey")

tsdens(cbind(ages, d18Opts0.2), "blue")
axis(1)
axis(2)

#compare auto-correlation terms between the models with different sized time bins
plot(density(p_18Ocal$BUGSoutput$sims.list$t_s.eps.ac), xlim=c(0,0.4), ylim=c(0,30),
     xlab = "auto-correlation term", main ="", col = "red")
lines(density(p_18Ocal0.2$BUGSoutput$sims.list$t_s.eps.ac),col = "purple") 
legend(0.3,30,c("0.5 Ma","0.2 Ma"),lwd=c(1,1),col=c("red","purple"))
#smaller time bins, less auto-correlation, wider uncertainty bounds for temp estimates
