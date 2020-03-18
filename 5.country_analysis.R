###
# Countries to analyze for EL: Rwanda, Niger, India, S.Africa = RWA, NER, IND, ZAF
###


# Functions for country analysis (assume X0 and del are given)

GetScaler <- function(cty.data) {
  gini.base  = cty.data$gini.base/100  # WDI has pct values.
  avg.base   = cty.data$avg.base  
  min.base   = cty.data$min.base  
  dle.thres  = cty.data$dle.thres 
  
  scl = list()
  for (dgini in seq(0.02, 0.2, 0.02)) {
    scl[[as.character(gini.base-dgini)]] = 
      TransformDistr(gini.base, gini.base-dgini, avg.base, min.base, dle.thres, 10)
  }
  sc = sapply(scl, function(x) x$k)
  
  return(list(data=cty.data, sc=sc))
}

GetRefLognorm <- function(cty.data) {
  gini.base  = cty.data$data$gini.base/100  # WDI has pct values.
  avg.base   = cty.data$data$avg.base  
  min.base   = cty.data$data$min.base 
  
  # X domain dependent on country distribution
  del = avg.base*1e3/12000 # x interval 
  x0 = seq(0, avg.base*1e3, del)
  
  y0 = RefLognorm(x0, gini.base, avg.base, min.base) 
  
  cty.data$y0=y0
  
  return(cty.data)
}


GetDF <- function(cty.data) {
  sc = cty.data$sc
  y0 = cty.data$y0
  dle.thres = cty.data$data$dle.thres 
  min.base  = cty.data$data$min.base  
  avg.base  = cty.data$data$avg.base  
  
  del = avg.base*1e3/12000 # x interval 
  x0 = seq(0, avg.base*1e3, del)
  
  xmat = sapply(sc, function(s) x0*s)
  ymat = sapply(sc, function(s) y0/s)
  
  df = data.frame(x0, y0)  
  # Append other PDF vectors to the DF
  for (i in 1:length(sc)) {
    idx = 2+i
    df[,idx] = approx(xmat[,i], ymat[,i], x0)$y  
    names(df)[idx] = paste0("y", i)
  }
  
  df$y0 = ShiftVector(df$y0, min.base/del)   # y0 is the baseline distribution
  df[, -c(1, 2)] = apply(df %>% select(y1:last_col()), 2, ShiftVector, dle.thres/del)
  
  cty.data$df = df
  return(cty.data)
}


# Find shares around DLE threshold
GetBaselinePDF <- function(cty.data) { # For the baseline PDF and for the income, d=0 & sc=1
  
  gini.base  = cty.data$data$gini.base/100  # WDI has pct values.
  avg.base   = cty.data$data$avg.base  
  min.base   = cty.data$data$min.base  
  dle.thres  = cty.data$data$dle.thres 
  # sc         = cty.data$sc  
  x0         = cty.data$x0
  
  # # Baseline distributions to consider for integration (!= reference when min.base>0)
  # f1 <- function(x, gini.base, avg.base, min.base, d=0, sc=1) { 
  #   d0 = min.base*sc + d  # Shift from min.ref
  #   x0 <- (x-d0)/sc 
  #   d1 <- RefLognorm(x0, gini.base, avg.base, min.base)/sc
  #   return(d1) # returns the prob density at x
  # }
  # xf1 <- function(x, gini.base, avg.base, min.base, d=0, sc=1) x*f1(x, gini.base, avg.base, min.base, d, sc)
  # dxf1 <- function(x, gini.base, avg.base, min.base, d=0, sc=1, offset) (offset-x)*f1(x, gini.base, avg.base, min.base, d, sc)
  
  # Share of population under DLE
  share.subDLE <- integrate(f1, min.base, dle.thres, 
                            gini.base, avg.base, min.base) # Base year (d=0, sc=1), share of sub-DLE population
  
  # Population avg = avg.base
  avg.calculated <- integrate(xf1, min.base, Inf, 
                              gini.base, avg.base, min.base) # Base year (d=0, sc=1), total wealth of sub-DLE population (normalized)
  
  # SubDLE avg
  wealth.subDLE <- integrate(xf1, min.base, dle.thres, 
                             gini.base, avg.base, min.base) # Base year (d=0, sc=1), total wealth of sub-DLE population (normalized)
  mean.subDLE <- wealth.subDLE$value / share.subDLE$value             # Base year, mean per capita of sub-DLE population
  
  # SubDLE gap
  gap.subDLE <- integrate(dxf1, min.base, dle.thres, 
                          gini.base, avg.base, min.base, offset=dle.thres) # Base year (d=0, sc=1), total wealth of sub-DLE population (normalized)
  mean.gap.subDLE <- gap.subDLE$value / share.subDLE$value             # Base year, mean per capita of sub-DLE population
  
  result = list(share.subDLE = share.subDLE, 
                avg.calculate = avg.calculated, 
                wealth.subDLE = wealth.subDLE, 
                mean.subDLE = mean.subDLE, 
                gap.subDLE = gap.subDLE, 
                mean.gap.subDLE = mean.gap.subDLE)
  
  return(result)
  # Find top x% which can offset the sub-DLE requirement 
  # What do we find? (For now, manually found)
  # (wealth.subDLE = wealth.rich) ->  No!
  # (gap.subDLE = surplus.rich)   ->  Yes!
  
  # thres.rich <- 80.76    # gap.subDLE = surplus.rich
  # # 104.6  # wealth.subDLE = wealth.rich
  # wealth.rich <- integrate(xf1, thres.rich, Inf, 
  #                          gini.base, avg.base, min.base)   # Manually find x=104.6, which can compensate all sub-DLE pop with the DLE threshold value
  # surplus.rich <- integrate(dxf1, thres.rich, Inf, 
  #                           gini.base, avg.base, min.base, offset=thres.rich)   # Manually find x=104.6, which can compensate all sub-DLE pop with the DLE threshold value
  # share.rich <- integrate(f1, thres.rich, Inf, 
  #                         gini.base, avg.base, min.base) 
  # mean.rich <- wealth.rich$value / share.rich$value 
  
}


# Plot country baselines
PlotCountry <- function(cty.data) {
  gini.base  = cty.data$data$gini.base/100  # WDI has pct values.
  avg.base   = cty.data$data$avg.base  
  min.base   = cty.data$data$min.base  
  dle.thres  = cty.data$data$dle.thres 
  
  del = avg.base*10/12000 # x interval 
  x0 = seq(0, avg.base*10, del)
  
  par(mar=c(2,2,2,2))
  # par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(x= seq(0, avg.base*10, del), y=RefLognorm(x0, gini.base, avg.base, min.base))
  title(cty.data$data$iso3c)
}


####
master.sub <- master %>% filter(iso3c %in% c('RWA', 'NER', 'IND', 'ZAF')) %>%
  group_by(iso3c) %>% 
  summarise_all(first) # Keep only the most recent obs.

master.sub <- master.sub %>% rename(gini.base = Gini, avg.base = GDP.PCAP) %>%
  mutate(min.base = 0, dle.thres = 1.9*365) %>% select(-(sd:ratio), -country, -region)

country.list <- split(master.sub, seq(nrow(master.sub)))
names(country.list) <- master.sub$iso3c

sc.list <- lapply(country.list, GetScaler)
y0.list <- lapply(sc.list, GetRefLognorm)
df.list <- lapply(y0.list, GetDF)
result.list <- lapply(df.list, GetBaselinePDF)

# Plot country baselines
lapply(df.list, PlotCountry)



# Progress 3/18/2020
# The top-to-bottom X vector decile cuts are resolved by adding jitters.
# Now the 4-country data are being used to generate baseline distributions.
# Need to find (manually for now) the 'rich' thresholds for each country.
# Also derive and compare gini changes for baseline, affine, top-to-bottom cases.
# Which affine (which sc value) to use? with more growth? (largest sc) <- sounds right
