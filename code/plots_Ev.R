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

d$ai = ceiling((70 - d$age) / ages.bin)

#compared to GMST stack by Westerhold
load("out/postTemp.rda") #make sure to run "driver_T.R" first

tp = p$BUGSoutput$sims.list$t_m
tp = tp[,-(ncol(tp))]
ages = ages[-length(ages)]#remove the last one -0.25 from the vector

tpts = apply(tp, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
tpts = t(tpts)

plot(-10, 0, xlab="Age (Ma)", ylab = "GMST (relative to preindustrial)", 
     xlim=c(70,0), ylim = c(-5,25), axes = FALSE)

points(d$age, d$temp, cex=0.5, col = "dark grey")

tsdens(cbind(ages, tpts), "firebrick4")
axis(1)
axis(2)