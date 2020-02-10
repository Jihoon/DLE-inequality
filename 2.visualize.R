require(reldist)
require(ineq)
require(stats)
require(ggplot2)
require(tidyverse)
require(EnvStats)
library(stats)

# mean_ene = 20 
# mlog=seq(0, log(mean_ene), log(mean_ene)/10) # For various means
# sdlog = sqrt(2*(log(mean_ene) - mlog))

# India DLE case

# Parameters are taken from 'Sample calculation - gini shift.xlsx' at \IIASA\DLE - Documents\WS3 - Documents\DLE and inequality
given_gini = 0.55 # Current Gini (energy)
avg_ene = 18.1    # Current avg GJ/cap
dle.thres = 11.3  # DLE GJ/cap
actual.min = 4.63  # Current minimum GJ/cap  

distr = "lognormal" # "pareto" # 

sc0 = 0.7
sc1 = 1.307
sc2 = 1.511
sc3 = 1.757
sc4 = 2.279

# Interval for the x values
del = 0.01
xval0 = seq(0, 120, del)
# xval0_org = xval0

if (distr == "lognormal") {
  gini_baselognormal = given_gini *avg_ene / (avg_ene-actual.min) # Assume scaler=1 for this base transformation

  sdlog_hat = qnorm((gini_baselognormal+1)/2)*sqrt(2) 
  mlog_hat = log(avg_ene-actual.min) - (sdlog_hat^2)/2 
  # sdlog_hat = qnorm((given_gini+1)/2)*sqrt(2) 
  # mlog_hat = log(avg_ene) - (sdlog_hat^2)/2 
  
  # Construct a base lognormal distribution for log x0 ~ N(mlog_hat, sdlog_hat)
  x0 = dlnorm(xval0, mlog_hat, sdlog_hat)
} else if (distr == "pareto") {
  # Pareto
  # dd = 10
  # x0 = dpareto(xval0, actual.min, shape=((1/given_gini)+1)/2)   # Gini = 1/(2*alpha-1)= 0.54 Staus quo
  source("3.Pareto.R")
}


# 
# x1 = dpareto(xval0, actual.min, shape=((1/given_gini)+1)/2/2)   # Gini = 1/(2*alpha-1)
# x2 = dpareto(xval0, actual.min*2, shape=((1/given_gini)+1)/2)   # Gini = 1/(2*alpha-1)
# xd = dpareto(xval0, actual.min+dd, shape=((1/given_gini)+1)/2)   # Gini = 1/(2*alpha-1)

# gini(2*r0)
# plot(xval0, x0, xlim=c(0,50), ylim=c(0, 0.5))
# plot(xval0, x1, xlim=c(0,50), ylim=c(0, 0.5))
# plot(xval0, x2, xlim=c(0,50), ylim=c(0, 0.5))
# plot(xval0, x0*2^(((1/given_gini)+1)/2), xlim=c(0,50), ylim=c(0, 0.5))
# plot(xval0*2, x0/2, xlim=c(0,50), ylim=c(0, 0.5))
# plot(xval0+dd, x0, xlim=c(0,50), ylim=c(0, 0.5))
# plot(xval0, xd, xlim=c(0,50), ylim=c(0, 0.5))

# # Test gini for pareto
# r0 = rpareto(1e7, actual.min, shape=((1/given_gini)+1)/2)   # Gini = 1/(2*alpha-1)
# r1 = r0*sc1 + dle.thres - sc1*actual.min
# r2 = r0*sc2 + dle.thres - sc2*actual.min
# r3 = r0*sc3 + dle.thres - sc3*actual.min
# r4 = r0*sc4 + dle.thres - sc4*actual.min
# r5 = r0*sc0 + dle.thres - sc0*actual.min
# print(paste(gini(r0), gini(r1), gini(r2), gini(r3), gini(r4), gini(r5)))


# Base scaled (x-axis)
if (distr == "pareto") {
  # xval0s = xval0
  xval0s = xval0*k_sol+d_sol
}
xval1 = xval0*sc1
xval2 = xval0*sc2
xval3 = xval0*sc3
xval4 = xval0*sc4
xval5 = xval0*sc0

# Base density scaled (Y-axis)
if (distr == "pareto") {
  # x0s = x0
  x0s = x0/k_sol
}
x1 = x0/sc1
x2 = x0/sc2
x3 = x0/sc3
x4 = x0/sc4
x5 = x0/sc0

# Create base DF with two base vectors


# Harmonize to the same X domain
if (distr == "pareto") {
  # xval0s = xval0
  df = data.frame(xval0, x0 = approx(xval0s, x0s, xval0)$y)
} else if (distr == "lognormal") {
  df = data.frame(xval0, x0)  
}
df$x1 = approx(xval1, x1, xval0)$y
df$x2 = approx(xval2, x2, xval0)$y
df$x3 = approx(xval3, x3, xval0)$y
df$x4 = approx(xval4, x4, xval0)$y
df$x5 = approx(xval5, x5, xval0)$y

# Introduce the shift (X-axis) - only for log-normal

if (distr == "lognormal") {
  df$x0 = c(rep(NA, actual.min/del), df$x0[1:(length(xval0)-actual.min/del)]) # A baseline distribution (lognormal x0 now with a shift)
  df$x1 = c(rep(NA, dle.thres/del),  df$x1[1:(length(xval0)-dle.thres/del)])
  df$x2 = c(rep(NA, dle.thres/del),  df$x2[1:(length(xval0)-dle.thres/del)])
  df$x3 = c(rep(NA, dle.thres/del),  df$x3[1:(length(xval0)-dle.thres/del)])
  df$x4 = c(rep(NA, dle.thres/del),  df$x4[1:(length(xval0)-dle.thres/del)])
  df$x5 = c(rep(NA, dle.thres/del),  df$x5[1:(length(xval0)-dle.thres/del)])
} else if (distr == "pareto") {
  # For pareto, even x0 needs an affine transformation.
  shift1 = as.integer((dle.thres - sc1*actual.min)/del)
  shift2 = as.integer((dle.thres - sc2*actual.min)/del)
  shift3 = as.integer((dle.thres - sc3*actual.min)/del)
  shift4 = as.integer((dle.thres - sc4*actual.min)/del)
  shift5 = as.integer((dle.thres - sc0*actual.min)/del)
  df$x1 = c(rep(NA, shift1), df$x1[1:(length(xval0)-shift1)])
  df$x2 = c(rep(NA, shift2), df$x2[1:(length(xval0)-shift2)])
  df$x3 = c(rep(NA, shift3), df$x3[1:(length(xval0)-shift3)])
  df$x4 = c(rep(NA, shift4), df$x4[1:(length(xval0)-shift4)])
  df$x5 = c(rep(NA, shift5), df$x5[1:(length(xval0)-shift5)])
}


gtext = data.frame(xv = apply(df[-1], 2, function(x){xval0[which(x==max(x, na.rm=T))]}),
                   yv = apply(df[-1], 2, function(x) max(x, na.rm=T)+0.001),
                   c = paste0('x', 0:5),
                   lab = c(paste0("G[x]:", given_gini),
                           paste0("list(G[z]:0.45, r:5.8*'%'/yr)"),
                           paste0("list(G[z]:0.48, r:6.8*'%'/yr)"),
                           paste0("list(G[z]:0.50, r:7.8*'%'/yr)"),
                           paste0("list(G[z]:0.54, r:9.8*'%'/yr)"),
                           paste0("list(G[z]:0.33, r:2.4*'%'/yr)")))
pdf(file = paste0("India GINI illustration-", distr, ".pdf"), width = 10, height = 6)
# Plot for Princeton w/s, May 2019
ggplot(df %>% select(-x2) %>% gather(key=dist, value=f, -xval0), aes(x=xval0, y=f, color=dist)) +
  geom_line(size=1) +
  geom_vline(xintercept = actual.min, linetype="dashed", colour="blue") +
  geom_vline(xintercept = dle.thres, linetype="dashed", colour="brown") +
  geom_text(aes(x=actual.min, label=paste0("Actual minimum : ", actual.min, "\n"), y=0.02), colour="blue", hjust=0.0, angle=90, size=4.5) +
  geom_text(aes(x=dle.thres, label=paste0("Desired minimum : ", dle.thres, "\n"), y=0.02), colour="brown", hjust=0.0, angle=90, size=4.5) +
  labs(x="Final energy consumption (GJ/cap)", y="Density")+
  # geom_text(aes(x=xval0[which(x0==max(x0, na.rm=T))], label=paste0("g_x: ", given_gini),   y=max(x0, na.rm=T)+0.001, colour="x0"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=xval0[which(x1==max(x1, na.rm=T))], label=paste0("g_z: 0.45, r: 5.8%"),  y=max(x1, na.rm=T)+0.001, colour="x1"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # # geom_text(aes(x=xval0[which(x2==max(x2, na.rm=T))], label=paste0("g: 0.48, r: 6.8%"),  y=max(x2, na.rm=T)+0.001, colour="x2"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=xval0[which(x3==max(x3, na.rm=T))], label=paste0("g_z: 0.50, r: 7.8%"),  y=max(x3, na.rm=T)+0.001, colour="x3"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=xval0[which(x4==max(x4, na.rm=T))], label=paste0("g_z: 0.54, r: 9.8%"), y=max(x4, na.rm=T)+0.001, colour="x4"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=xval0[which(x5==max(x5, na.rm=T))], label=paste0("g_z: 0.33, r: 2.4%")), y=max(x5, na.rm=T)+0.001, colour="x5"), size=4, hjust=0.0, vjust=0.0, data=df) +
  geom_text(data=gtext[-3,], aes(x=xv, label=lab, y=yv, colour=c), size=4, hjust=0.0, vjust=0.0, parse=TRUE) +
  theme(legend.position = "none") +
  xlim(c(0, 30))

dev.off()



### Try to integrate

# Base lognormal
f0 <- function(x) dlnorm(x, mlog_hat, sdlog_hat)
fp0 <- function(x) dpareto(x, 5, 1) # test
g0 <- function(x) x*f0(x)
gp0 <- function(x) x*fp0(x)

# Focus on top 5% population
topinc <- 0.05

# find the energy threshold for the high income population (at the base lognormal distr)
xval <- qlnorm(1-topinc, mlog_hat, sdlog_hat) # check the top x%

# test integrate function
integrate(f0, 0, xval)
integrate(g0, 0, Inf) # mean
integrate(g0, xval, Inf)$value/integrate(f0, xval, Inf)$value # mean of the top 5%

# Current IND distr with shift = actual.min > 0 (Gini=0.55 (given_gini))
f1 <- function(x, d=actual.min, sc=1) { # sc=1 assumed 
  x0 <- (x-d)/sc
  d1 <- dlnorm(x0, mlog_hat, sdlog_hat)/sc
  return(d1) # returns the prob density 
}
g1 <- function(x, d=actual.min, sc=1) x*f1(x, d, sc)  

f2 <- function(x, d, sc) { # d here means the shift from f1 distribution (dle.thres = sc*actual.min + d)
  x0 <- (x-d)/sc
  d1 <- f1(x0)/sc
  return(d1)
}
g2 <- function(x, d, sc) x*f2(x, d, sc)  

integrate(f1, 0, xval+actual.min, d=actual.min, sc=1)

# Current avg wealth of upper 10% 
top_mean_now <- integrate(g1, xval+actual.min, Inf, d=actual.min, sc=1)$value/integrate(f1, xval+actual.min, Inf, d=actual.min, sc=1)$value

# integrate(f2, 0, Inf, d=dle.thres-sc1*actual.min, sc=sc1) # gives 1

# integrate(f2, xval*sc1+dle.thres, Inf, d=dle.thres-sc1*actual.min, sc=sc1)
integrate(g2, xval*sc1+dle.thres, Inf, d=dle.thres-sc1*actual.min, sc=sc1)$value/integrate(f2, xval*sc1+dle.thres, Inf, d=dle.thres-sc1*actual.min, sc=sc1)$value
# integrate(f2, xval*sc2+dle.thres, Inf, d=dle.thres-sc2*actual.min, sc=sc2)
integrate(g2, xval*sc2+dle.thres, Inf, d=dle.thres-sc3*actual.min, sc=sc2)$value/integrate(f2, xval*sc2+dle.thres, Inf, d=dle.thres-sc2*actual.min, sc=sc2)$value
# integrate(f2, xval*sc3+dle.thres, Inf, d=dle.thres-sc3*actual.min, sc=sc3)
integrate(g2, xval*sc3+dle.thres, Inf, d=dle.thres-sc3*actual.min, sc=sc3)$value/integrate(f2, xval*sc3+dle.thres, Inf, d=dle.thres-sc3*actual.min, sc=sc3)$value

# integrate(f2, xval*sc0+dle.thres, Inf, d=dle.thres-sc0*actual.min, sc=sc0)
sc0 = 0.7 # we can test different sc values.
top_mean_dle <- integrate(g2, xval*sc0+dle.thres, Inf, d=dle.thres-sc0*actual.min, sc=sc0)$value/integrate(f2, xval*sc0+dle.thres, Inf, d=dle.thres-sc0*actual.min, sc=sc0)$value
c(top_mean_now, top_mean_dle, sc0)

###

x1=rlnorm(1e6, mlog_hat, sdlog_hat)
a1 <- Lc(x1, plot=FALSE)


x2=x1*sc1+(dle.thres-sc1*actual.min)
a2 <- Lc(x2, plot=FALSE)
Gini(x2)

x3=x1*sc2+(dle.thres-sc2*actual.min)
a3 <- Lc(x3, plot=FALSE)
Gini(x3)

x4=x1*sc3+(dle.thres-sc3*actual.min)
a4 <- Lc(x4, plot=FALSE)
Gini(x4)

x5=x1*sc4+(dle.thres-sc4*actual.min)
a5 <- Lc(x5, plot=FALSE)
Gini(x5)

df = rbind(data.frame(n=1:length(x1), dist="base", x=x1),
           data.frame(n=1:length(x1), dist="case1", x=x2),
           data.frame(n=1:length(x1), dist="case2", x=x3),
           data.frame(n=1:length(x1), dist="case3", x=x4),
           data.frame(n=1:length(x1), dist="case4", x=x5))

plot( density(x1, bw="nrd"), col='red', ylim=c(0,0.08), xlim=c(0, 150), main="", lwd=2)
lines(density(x2), col='blue', lwd=2)
lines(density(x3), col='green', lwd=2)
lines(density(x4), col='purple', lwd=2)
lines(density(x5), col='brown', lwd=2)

hist(x1, breaks=100000, xlim=c(0, 150))

ggplot(df, aes(x=x, color=dist)) + 
  geom_density() +
  lims(x=c(0, 150))


# data.frame(density(x1), x2=density(x2), x3=density(x3), x4=density(x4))

plot(a1, col='red')
lines(a2, col='blue')
lines(a3, col='green')
lines(a4, col='purple')
lines(a5, col='brown')





#==============================

# Making sense of lnorm functions
x0 = rlnorm(50000,0,.6)
x1 = rlnorm(50000,0,.8)
x2 = rlnorm(50000,0,.1)
lapply(list(x0, x1, x2), summary)
grid = seq(0,25,.1)

plot(grid,dlnorm(grid, 0, .6),type="l",xlab="x",ylab="f(x)")  # PDF
lines(dlnorm(grid,0,1),col="red")
lines(dlnorm(grid,0,2),col="green")

plot(grid,dlnorm(grid, 0, .6),type="l",xlab="x",ylab="f(x)")  # PDF
lines(density(x0),col="red")
lines(density(x1),col="green")
lines(density(x2),col="blue")

plot(grid,plnorm(grid, 1, .6),type="l",xlab="x",ylab="F(x)")  # CDF

lines(density(x),col="red")

legend("topright",c("True Density","Estimate"),lty=1,col=1:2)


### Deriving Ginis for given lnorm mean and min (=0 for now).
# Keep same mean for lognormal above zero, then what can be possible ginis?
mean_ene = 20 
mlog=seq(0, log(mean_ene), log(mean_ene)/10) # For various means
sdlog = sqrt(2*(log(mean_ene) - mlog))
gini = 2*pnorm(sdlog/sqrt(2))-1

for (i in 1:length(mlog)) {
  # plot(dlnorm(grid, mlog[i], sdlog[i]),type="l",xlab="x",ylab="f(x)")  
  x=rlnorm(1e6, mlog[i], sdlog[i])
  # plot(density(x), type="l",xlab="x",ylab="f(x)", xlim=c(1, 50))  
  hist(x, breaks=10000, 
       xlim=c(0, mean_ene*5), #ylim=c(0, 1500), 
       freq=FALSE)
  abline(v=mean(x))
  title(sub=gini[i])
  print(summary(x))
  print(paste("Gini =",gini[i]))
  # par(new=TRUE)
}
dev.off()

# Compare above with the actual India Gini in 2011 (0.35)
# 0.53 - BRA, 0.63 - ZAF
# With given mean and gini
SpanEnergy <- function(gn, avg_ene) {
  given_gini = gn #0.35
  max_avg_ene = 40
  sdlog_hat = qnorm((given_gini+1)/2)*sqrt(2)  # 
  mlog_hat = log(avg_ene) - (sdlog_hat^2)/2 
  
  sdlog_new = sdlog_hat * log(max_avg_ene-avg_ene)
  c = 1
  mlog_new = c
  
  mlog_hat * log(max_avg_ene-avg_ene)
  
  # mlog_given1 = log(max_avg_ene) - (sdlog_given^2)/2 
  # + (40-avg_ene)
  x1=rlnorm(1e7, mlog_hat, sdlog_hat) 
  x2=rlnorm(1e7, c, sdlog_new) - (max_avg_ene-avg_ene)
  # print(c(summary(x)[1], summary(x)[6]))
  # print(gini(x))
  hist(x1, breaks=10000, 
       # xlim=c(0, avg_ene*5), 
       xlim=c(0, 200),#ylim=c(0, 1500), 
       freq=FALSE)
  abline(v=mean(x1))
  print(summary(x1))
  hist(x2, breaks=10000, 
       # xlim=c(0, avg_ene*5), 
       xlim=c(0, 200),#ylim=c(0, 1500), 
       freq=FALSE)
  abline(v=mean(x2))  
  print(summary(x2))
  print(paste("Gini1 =", Gini(x1)))
  print(paste("Gini2 =", Gini(x2)))
  
  return(list(x1, x2))
}

list[x1, x2] <- SpanEnergy(0.5, 35)
mean_vec = seq(5, 40, 5)
par(mfrow = c(4, 2))
x_income <- sapply(mean_vec, SpanEnergy, gn=0.5)

# SpanEnergy(0.35, 3)


SpanEnergy2(0.35, 20)


x=(rlnorm(1e7, 1, 0.5))
mean(x)
Gini(x)

SpanEnergy3 <- function(xmin, avg_ene) { #avg_ene : real avg (with shift)
  
  mean.base = avg_ene - xmin # Assume we want to achieve min=35, avg_ene=40 (gap=35) -> mean.base
  # gini.goal = gn
  # sd.org = qnorm((gini.goal+1)/2)*sqrt(2) # SD When there is no shift.
  
  # Get mu and sd giving the same mean.base
  mu.range = log(mean.base) + c(-20)#, -10, -7, -5, -3, -1, -0.5, -0.3, -0.1, -0.01)
  sd.range = sqrt(2*(log(mean.base) - mu.range))
  
  for (i in 1:length(mu.range)) {
    x = (rlnorm(5e7, mu.range[i], sd.range[i]))  
    # plot(x)
    print(paste("Gini", Gini(x), Gini(x+xmin), 2*pnorm(sd.range[i]/sqrt(2))-1))
    print(paste("Mean", mean(x), mean(x+xmin), exp((sd.range[i])^2/2+mu.range[i])))
    print(paste("mu/sd", mu.range[i], sd.range[i]))
  }
}

SpanEnergy3(8, 15)



x=(rlnorm(1e7, 2, 0.5)) 
hist(x, breaks=10000, freq=FALSE)

# plot(dlnorm(seq(0,100000,.01), 0, 0.01))
# x = rlnorm(500,1,.6)
grid = seq(0,100000,.01)
s = 1000
plot(grid,dlnorm(grid, 0, s),type="l",xlab="x",ylab="f(x)")
plot(grid,plnorm(grid, 0, s),type="l",xlab="x",ylab="f(x)")

# Horizontal shift
lx = log(x+8) 
hist(lx, breaks=10000, freq=FALSE) # This is not normal.
Gini(x+8)
Gini(x)


lx.shifted = lx*log(7) #log(x+7)
hist(lx.shifted, breaks=10000, freq=FALSE)

x.sh=(rnorm(1e7, 2*log(7), 0.5*log(7))) 
hist(x.sh, breaks=10000, freq=FALSE)
# x=(rlnorm(1e7, 200, 0.5)) 

Gini(x)


hist(x+5, breaks=10000, xlim=c(0,max(x)+5),
     freq=FALSE)

# Testing pareto distribution
# library(rmutil) # or
library(EnvStats)
xmin = 8
xmean = 14
alpha = xmean/(xmean-xmin)
gini = 1/(2*alpha -1)

xmin = 15
v.mean = seq(xmin+1, xmin+50, 1)
v.alpha = v.mean/(v.mean-xmin)
v.gini = 1/(2*v.alpha -1)
v.pareto = data.frame(v.mean, v.gini, v.alpha)

ggplot(v.pareto, aes(x=v.mean, y=v.gini)) +
  geom_point()

sprintf("alpha=%f, gini=%f", alpha, gini)
# nsample
x.p <- ppareto(1:10000, xmin, alpha)
d.p <- dpareto(1:10000, xmin, alpha)
r.p <- rpareto(100000, xmin, alpha)


Gini(r.p)
mean(r.p)
alpha*xmin/(alpha-1)
hist(r.p, breaks="FD", freq=T, xlim=c(xmin,200))
plot(x.p)
plot(d.p)
(d.p)[1:100]


# I derived a more general formula for gini calculation after the shift/scaling.
# This is for testing the furmula 
# Just shift the given gini (before shift) lognormal to fixed mean after shift
SpanEnergy4 <- function(gn, avg_ene) { #avg_ene : real avg (with shift)
  # gaps = seq(0, avg_ene-5, 5)
  
  # for (gap in gaps) {
  # mean.base = avg_ene - gap # Assume we want to achieve min=35, avg_ene=40 (gap=35) -> mean.base
  # gini.goal = gn
  
  shifts = setdiff(seq(0, 1, 0.25), c(0,1))
  par(mfrow = c(length(shifts)+1, 1))
  
  sd.org = qnorm((gn+1)/2)*sqrt(2) # SD When there is no shift.
  
  mu = log(avg_ene) - (sd.org^2)/2 
  # if (gap==0) {
  x.org = (rlnorm(1e5, mu, sd.org))
  hist(x.org, breaks=10000, freq=FALSE, xlim=c(0, 100))
  abline(v=mean(x.org))  
  title(sub=paste('Gini=', Gini(x.org)))
  print(summary(x.org))
  #   next
  # }    # sd.min = sqrt(log(mean.base)*2)
  
  for (shift in shifts) {
    x.new = x.org*(1-shift) + avg_ene*shift
    print(paste('Gini.new =', Gini(x.new), 'Gini.old =', Gini(x.org),
                '; Mean.new =', mean(x.new), '; shift =', shift))
    
    hist(x.new, breaks=1000, freq=FALSE, xlim=c(0, 100))
    abline(v=mean(x.new))  
    title(sub=paste('Gini=', Gini(x.new)))
    print(summary(x.new))
  }
  # }
}

SpanEnergy4(0.4, 40)
# For any distributions
SpanEnergy5 <- function(avg_ene, dist="lognormal") { #avg_ene : real avg (with shift), dist="beta" or "pareto"
  
  scales = setdiff(seq(0, 1, 0.2), c(0,1))
  par(mfrow = c(length(scales)+1, 1), mar=c(1,1,1,1))
  
  n=10000
  if (dist=="beta") {
    x.org = rbeta(n, 2, 5)
  } else if (dist=="pareto") {
    x.org = rpareto(n, 1, 1) # shape (alpha) = 1
  } else {
    x.org = rlnorm(n)
  }
  
  hist(x.org, breaks=10000, freq=FALSE, xlim=c(0, as.numeric(quantile(x.org, 0.85))))
  abline(v=mean(x.org))  
  title(sub=paste('Gini=', Gini(x.org)))
  print(paste(dist, summary(x.org)))
  
  for (scaler in scales) {
    x.new = x.org*scaler + (avg_ene-mean(x.org)*scaler)
    print(paste('Gini.new =', Gini(x.new), 'Gini.old =', Gini(x.org),
                '; Mean.new =', mean(x.new), '; scaler =', scaler, 
                '; scaler (der) =', scaler*mean(x.org)/(scaler*mean(x.org) + (avg_ene-scaler*mean(x.org))),
                '; scaler (calc) =',  Gini(x.new)/Gini(x.org)))
    
    hist(x.new, breaks=10000, freq=FALSE, xlim=c(0, as.numeric(quantile(x.new, 0.85))))
    abline(v=mean(x.new))  
    title(sub=paste('Gini=', Gini(x.new)))
    print(summary(x.new))
  }
}


SpanEnergy5(4, "beta")
SpanEnergy5(40, "pareto")
SpanEnergy5(40, "lognormal")
dev.off()

