require(reldist)
require(ineq)
require(stats)
require(ggplot2)
require(tidyverse)
require(EnvStats)
library(stats)
library(devtools)
source_url("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")


# India DLE case

# Parameters are taken from 'Sample calculation - gini shift.xlsx' at \IIASA\DLE - Documents\WS3 - Documents\DLE and inequality
gini.base = 0.55 # Current Gini (energy)
avg.base = 18.1    # Current avg GJ/cap
dle.thres = 11.3  # DLE GJ/cap
min.base = 4.6236  # Current minimum GJ/cap  

distr = "lognormal" # "pareto" # 

# sc = c(0.7,
#        1.307,
#        1.511,
#        1.757,
#        2.279)


# Derive scaling factors (k) for the target (delta) Gini values
scl = list()
for (dgini in seq(0, 0.2, 0.02)) {
  scl[[as.character(gini.base-dgini)]] = TransformDistr(gini.base, gini.base-dgini, avg.base, min.base, dle.thres, 10)
}
sc = sapply(scl, function(x) x$k)


# Interval for the x domain
del = 0.01 # x interval 
x0 = seq(0, 120, del)


# Get the reference PDF values for the x domain (x0)
if (distr == "lognormal") {
  # list[y0, F0] = RefLognorm(x0, gini.base, avg.base, min.base) 
  y0 = RefLognorm(x0, gini.base, avg.base, min.base) 
} else if (distr == "pareto") {
  y0 = RefPareto(x0, gini.base, avg.base, min.base)  
}


# Scale the PDFs with the scalers (both for x and y axes)
###! Still requires update for Pareto
# X-axis
if (distr == "pareto") {
  x0s = x0*k_sol+d_sol
}
xmat = sapply(sc, function(s) x0*s)

# Y-axis
if (distr == "pareto") {
  y0s = y0/k_sol
}
ymat = sapply(sc, function(s) y0/s)


# Create base DF with two base vectors (x0, y0)
# Harmonize to the same X domain
###! Still requires update for Pareto
if (distr == "pareto") {
  df = data.frame(x0, y0 = approx(x0s, y0s, x0)$y)
} else if (distr == "lognormal") {
  df = data.frame(x0, y0)  
}


# Append other PDF vectors to the DF
for (i in 1:length(sc)) {
  idx = 2+i
  df[,idx] = approx(xmat[,i], ymat[,i], x0)$y  
  names(df)[idx] = paste0("y", i)
}


# Function to introduce the shift (X-axis) - only for log-normal
ShiftVector <- function(v, shift) {
  shift = as.integer(shift)
  xshift = c(rep(NA, shift),  v[1:(length(v)-shift)])
  return(xshift)
}


# Bring the shift to DF
if (distr == "lognormal") {
  df$y0 = ShiftVector(df$y0, min.base/del)   # y0 is the baseline distribution
  df[, -c(1, 2)] = apply(df %>% select(y1:last_col()), 2, ShiftVector, dle.thres/del)
  
} else if (distr == "pareto") {
  # For pareto, even x0 needs an affine transformation.
  for (i in 1:length(sc)) {
    sh = as.integer((dle.thres - sc[i]*min.base)/del)
    df[,i+2] = ShiftVector(df[,i+2], sh)
  }
}


gtext = data.frame(xv = apply(df[-1], 2, function(x){x0[which(x==max(x, na.rm=T))]}),
                   yv = apply(df[-1], 2, function(x) max(x, na.rm=T)+0.001),
                   c = paste0('x', 0:5),
                   lab = c(paste0("G[x]:", gini.base),
                           paste0("list(G[z]:0.45, r:5.8*'%'/yr)"),
                           paste0("list(G[z]:0.48, r:6.8*'%'/yr)"),
                           paste0("list(G[z]:0.50, r:7.8*'%'/yr)"),
                           paste0("list(G[z]:0.54, r:9.8*'%'/yr)"),
                           paste0("list(G[z]:0.33, r:2.4*'%'/yr)")))
pdf(file = paste0("India GINI illustration-", distr, ".pdf"), width = 10, height = 6)
# Plot for Princeton w/s, May 2019
ggplot(df %>% select(-y2) %>% gather(key=dist, value=f, -x0), aes(x=x0, y=f, color=dist)) +
  geom_line(size=1) +
  geom_vline(xintercept = min.base, linetype="dashed", colour="blue") +
  geom_vline(xintercept = dle.thres, linetype="dashed", colour="brown") +
  geom_text(aes(x=min.base, label=paste0("Actual minimum : ", min.base, "\n"), y=0.02), colour="blue", hjust=0.0, angle=90, size=4.5) +
  geom_text(aes(x=dle.thres, label=paste0("Desired minimum : ", dle.thres, "\n"), y=0.02), colour="brown", hjust=0.0, angle=90, size=4.5) +
  labs(x="Final energy consumption (GJ/cap)", y="Density")+
  # geom_text(aes(x=x0[which(x0==max(x0, na.rm=T))], label=paste0("g_x: ", gini.base),   y=max(x0, na.rm=T)+0.001, colour="x0"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x1==max(x1, na.rm=T))], label=paste0("g_z: 0.45, r: 5.8%"),  y=max(x1, na.rm=T)+0.001, colour="x1"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # # geom_text(aes(x=x0[which(x2==max(x2, na.rm=T))], label=paste0("g: 0.48, r: 6.8%"),  y=max(x2, na.rm=T)+0.001, colour="x2"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x3==max(x3, na.rm=T))], label=paste0("g_z: 0.50, r: 7.8%"),  y=max(x3, na.rm=T)+0.001, colour="x3"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x4==max(x4, na.rm=T))], label=paste0("g_z: 0.54, r: 9.8%"), y=max(x4, na.rm=T)+0.001, colour="x4"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x5==max(x5, na.rm=T))], label=paste0("g_z: 0.33, r: 2.4%")), y=max(x5, na.rm=T)+0.001, colour="x5"), size=4, hjust=0.0, vjust=0.0, data=df) +
  geom_text(data=gtext[-3,], aes(x=xv, label=lab, y=yv, colour=c), size=4, hjust=0.0, vjust=0.0, parse=TRUE) +
  theme(legend.position = "none") +
  xlim(c(0, 30))

dev.off()


### Integrate distributions from above

# The reference lognormal distr (not relevent for the integration analysis)
f0 <- RefLognorm
xf0 <- function(x, gini, avg, min) x*RefLognorm(x, gini, avg, min)


# Baseline distributions to consider for integration (!= reference when min.base>0)
f1 <- function(x, gini.base, avg.base, min.base, d, sc) { # Assume the target d & sc are derived from TransformDistr above
  d0 = min.base*sc + d  # Shift from min.ref
  x0 <- (x-d0)/sc 
  d1 <- RefLognorm(x0, gini.base, avg.base, min.base)/sc
  return(d1) # returns the prob density at x
}
xf1 <- function(x, gini.base, avg.base, min.base, d, sc) x*f1(x, gini.base, avg.base, min.base, d, sc)


# Attempt to define a truncated PDF with a dirac-delta at dle.thres level, 
# but not working (don't know why yet). So not used for now
f1.trunc <- function(x, f, min, max, w) { #, gini.base, avg.base, min.base, d, sc) { # Assume the target d & sc are derived from TransformDistr above
  
  library(distr)
  x = as.matrix(x)
  
  y = f1(x, gini.base=gini.base, avg.base=avg.base, min.base=min.base, d=0, sc=1)
  y = as.matrix(y)
  
  y[x>max] = 0
  y[x<min] = 0
  
  D = Dirac(x[which(y>0)[1]]) # Attempt to add Dirac, but not working
  y = y + d(D)(x)*w

  return(as.numeric(y)) # returns the prob density 
}



# a <- scl[[6]]

# Share of population under DLE
share.subDLE <- integrate(f1, min.base, dle.thres, 
                          gini.base, avg.base, min.base, d=0, sc=1) # Base year (d=0, sc=1), share of sub-DLE population

# Population avg = avg.base
avg.calculated <- integrate(xf1, min.base, Inf, 
                           gini.base, avg.base, min.base, d=0, sc=1) # Base year (d=0, sc=1), total wealth of sub-DLE population (normalized)

# SubDLE avg
wealth.subDLE <- integrate(xf1, min.base, dle.thres, 
                          gini.base, avg.base, min.base, d=0, sc=1) # Base year (d=0, sc=1), total wealth of sub-DLE population (normalized)
mean.subDLE <- wealth.subDLE$value / share.subDLE$value             # Base year, mean per capita of sub-DLE population


# Find top x% which can offset the sub-DLE requirement (wealth.subDLE = wealth.rich)
thres.rich <- 104.6  # For now, manually found
wealth.rich <- integrate(xf1, thres.rich, Inf, 
          gini.base, avg.base, min.base, d=0, sc=1)   # Manually find x=104.6, which can compensate all sub-DLE pop with the DLE threshold value
share.rich <- integrate(f1, thres.rich, Inf, 
          gini.base, avg.base, min.base, d=0, sc=1) 
mean.rich <- wealth.rich$value / share.rich$value 


# Draw the baseline distribution (to analyze top to bottom redistribution)
X <- DrawRefLognorm(1e5, gini.base, avg.base, min.base) + min.base
X[X<dle.thres] <- dle.thres
X[X>104.6] <- 104.6




# Test truncated distribution
share.mid <- integrate(f1, dle.thres, thres.rich, 
          gini.base, avg.base, min.base, d=0, sc=1)
w <- 1-share.mid$value

integrate(f1.trunc, 0, Inf, f1, min=dle.thres, max=thres.rich, w) # ERROR: the integral is probably divergent



a=f1.trunc(x0, f1, min=dle.thres, max=104.6, w)

D = Dirac(100)
d(D)(x0)[10001]
integrate(d(D), -Inf, Inf)




















# integrate(xf1, min.base, Inf, gini.base, avg.base, min.base, d=0, sc=1) # = min.base (test)
# integrate(xf1, 0, Inf, gini.base, avg.base, min.base, a$d, a$k) # = avg.new (test)
integrate(xf1, 0, Inf, gini.base, avg.base, min.base, a$d, a$k) # = avg.new

x.upper <- qlnorm(1-share.subDLE, mlog_hat, sdlog_hat) # check the top x%

integrate(f1, min.base, x.upper, gini.base, avg.base, min.base)






# Focus on top 5% population
topinc <- 0.05
x.upper <- qlnorm(1-topinc, mlog_hat, sdlog_hat) # check the top x%


#integrate(f0, 0, Inf, gini.base, avg.base, min.base) # =1
integrate(f0, min.base, x.upper, gini.base, avg.base, min.base)

### Try to integrate

# Base lognormal
f0 <- function(x) dlnorm(x, mlog_hat, sdlog_hat)
fp0 <- function(x) dpareto(x, 5, 1) # test
g0 <- function(x) x*f0(x)
gp0 <- function(x) x*fp0(x)

# find the energy threshold for the high income population (at the base lognormal distr)
xval <- qlnorm(1-topinc, mlog_hat, sdlog_hat) # check the top x%

# test integrate function
integrate(f0, 0, xval)
integrate(g0, 0, Inf) # mean
integrate(g0, xval, Inf)$value/integrate(f0, xval, Inf)$value # mean of the top 5%

# Current IND distr with shift = min.base > 0 (Gini=0.55 (gini.base))
f1 <- function(x, d=min.base, sc=1) { # sc=1 assumed 
  x0 <- (x-d)/sc
  d1 <- dlnorm(x0, mlog_hat, sdlog_hat)/sc
  return(d1) # returns the prob density 
}
g1 <- function(x, d=min.base, sc=1) x*f1(x, d, sc)  

f2 <- function(x, d, sc) { # d here means the shift from f1 distribution (dle.thres = sc*min.base + d)
  x0 <- (x-d)/sc
  d1 <- f1(x0)/sc
  return(d1)
}
g2 <- function(x, d, sc) x*f2(x, d, sc)  

integrate(f1, 0, xval+min.base, d=min.base, sc=1)

# Current avg wealth of upper 10% 
top_mean_now <- integrate(g1, xval+min.base, Inf, d=min.base, sc=1)$value/integrate(f1, xval+min.base, Inf, d=min.base, sc=1)$value

# integrate(f2, 0, Inf, d=dle.thres-sc1*min.base, sc=sc1) # gives 1

# integrate(f2, xval*sc1+dle.thres, Inf, d=dle.thres-sc1*min.base, sc=sc1)
integrate(g2, xval*sc1+dle.thres, Inf, d=dle.thres-sc1*min.base, sc=sc1)$value/integrate(f2, xval*sc1+dle.thres, Inf, d=dle.thres-sc1*min.base, sc=sc1)$value
# integrate(f2, xval*sc2+dle.thres, Inf, d=dle.thres-sc2*min.base, sc=sc2)
integrate(g2, xval*sc2+dle.thres, Inf, d=dle.thres-sc3*min.base, sc=sc2)$value/integrate(f2, xval*sc2+dle.thres, Inf, d=dle.thres-sc2*min.base, sc=sc2)$value
# integrate(f2, xval*sc3+dle.thres, Inf, d=dle.thres-sc3*min.base, sc=sc3)
integrate(g2, xval*sc3+dle.thres, Inf, d=dle.thres-sc3*min.base, sc=sc3)$value/integrate(f2, xval*sc3+dle.thres, Inf, d=dle.thres-sc3*min.base, sc=sc3)$value

# integrate(f2, xval*sc0+dle.thres, Inf, d=dle.thres-sc0*min.base, sc=sc0)
sc0 = 0.7 # we can test different sc values.
top_mean_dle <- integrate(g2, xval*sc0+dle.thres, Inf, d=dle.thres-sc0*min.base, sc=sc0)$value/integrate(f2, xval*sc0+dle.thres, Inf, d=dle.thres-sc0*min.base, sc=sc0)$value
c(top_mean_now, top_mean_dle, sc0)

###

x1=rlnorm(1e6, mlog_hat, sdlog_hat)
a1 <- Lc(x1, plot=FALSE)


x2=x1*sc1+(dle.thres-sc1*min.base)
a2 <- Lc(x2, plot=FALSE)
Gini(x2)

x3=x1*sc2+(dle.thres-sc2*min.base)
a3 <- Lc(x3, plot=FALSE)
Gini(x3)

x4=x1*sc3+(dle.thres-sc3*min.base)
a4 <- Lc(x4, plot=FALSE)
Gini(x4)

x5=x1*sc4+(dle.thres-sc4*min.base)
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
SpanEnergy <- function(gn, avg.base) {
  gini.base = gn #0.35
  max_avg.base = 40
  sdlog_hat = qnorm((gini.base+1)/2)*sqrt(2)  # 
  mlog_hat = log(avg.base) - (sdlog_hat^2)/2 
  
  sdlog_new = sdlog_hat * log(max_avg.base-avg.base)
  c = 1
  mlog_new = c
  
  mlog_hat * log(max_avg.base-avg.base)
  
  # mlog_given1 = log(max_avg.base) - (sdlog_given^2)/2 
  # + (40-avg.base)
  x1=rlnorm(1e7, mlog_hat, sdlog_hat) 
  x2=rlnorm(1e7, c, sdlog_new) - (max_avg.base-avg.base)
  # print(c(summary(x)[1], summary(x)[6]))
  # print(gini(x))
  hist(x1, breaks=10000, 
       # xlim=c(0, avg.base*5), 
       xlim=c(0, 200),#ylim=c(0, 1500), 
       freq=FALSE)
  abline(v=mean(x1))
  print(summary(x1))
  hist(x2, breaks=10000, 
       # xlim=c(0, avg.base*5), 
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

SpanEnergy3 <- function(xmin, avg.base) { #avg.base : real avg (with shift)
  
  mean.base = avg.base - xmin # Assume we want to achieve min=35, avg.base=40 (gap=35) -> mean.base
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
SpanEnergy4 <- function(gn, avg.base) { #avg.base : real avg (with shift)
  # gaps = seq(0, avg.base-5, 5)
  
  # for (gap in gaps) {
  # mean.base = avg.base - gap # Assume we want to achieve min=35, avg.base=40 (gap=35) -> mean.base
  # gini.goal = gn
  
  shifts = setdiff(seq(0, 1, 0.25), c(0,1))
  par(mfrow = c(length(shifts)+1, 1))
  
  sd.org = qnorm((gn+1)/2)*sqrt(2) # SD When there is no shift.
  
  mu = log(avg.base) - (sd.org^2)/2 
  # if (gap==0) {
  x.org = (rlnorm(1e5, mu, sd.org))
  hist(x.org, breaks=10000, freq=FALSE, xlim=c(0, 100))
  abline(v=mean(x.org))  
  title(sub=paste('Gini=', Gini(x.org)))
  print(summary(x.org))
  #   next
  # }    # sd.min = sqrt(log(mean.base)*2)
  
  for (shift in shifts) {
    x.new = x.org*(1-shift) + avg.base*shift
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
SpanEnergy5 <- function(avg.base, dist="lognormal") { #avg.base : real avg (with shift), dist="beta" or "pareto"
  
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
    x.new = x.org*scaler + (avg.base-mean(x.org)*scaler)
    print(paste('Gini.new =', Gini(x.new), 'Gini.old =', Gini(x.org),
                '; Mean.new =', mean(x.new), '; scaler =', scaler, 
                '; scaler (der) =', scaler*mean(x.org)/(scaler*mean(x.org) + (avg.base-scaler*mean(x.org))),
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

