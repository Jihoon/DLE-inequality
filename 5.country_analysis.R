###
# Countries to analyze for EL: Rwanda, Niger, India, S.Africa = RWA, NER, IND, ZAF
###

# Note: This file extends scripts in '2.distribution_analysis.R' into function forms,
#       so it can easily run the four-country calculations easily.


source("3.functions.R")

# Functions for country analysis (assume X0 and del are given)

GetScaler <- function(cty.data) {
  gini.base  = cty.data$gini.base/100  # WDI has pct values.
  avg.base   = cty.data$avg.base  
  min.base   = cty.data$min.base  
  dle.thres  = cty.data$dle.thres 
  
  scl = list()
  for (dgini in seq(0.01, 0.2, 0.01)) {
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
  dle.thres   = cty.data$data$dle.thres 
  
  # X domain dependent on country distribution
  del = dle.thres*1e3/12000 # x interval (based on dle.thres to make the shifted distribution start right from )
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
  
  # del should be defined based on dle.thres, to make the shifted distribution start exactly from dle.thres. (x0 coordinate) 
  del = dle.thres*1e3/12000 # x interval 
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
  
  
  # Get the rich share which can fill the subDLE gap
  eqn=function(th) {
    # This first integration is inverted (negative), so need to 'add'.
    # This integration frequently gives roundoff error with 'integrate', so use 'cubintegrate' instead.
    cubintegrate(f=dxf1, th, Inf, method = "pcubature", 
                 gini.b=gini.base, avg.b=avg.base, min.b=min.base, offset=th)$integral +
      integrate(dxf1, min.base, dle.thres, 
                gini.base, avg.base, min.base, offset=dle.thres)$value
    }

  # Solve the equation and find the solution for thres.rich
  # When dle.thres is too high, there will be countries with no solution.
  thres.rich = tryCatch(uniroot(eqn, c(avg.base, avg.base*1000))$root, 
                        error = function(e) {
                          print(paste(cty.data$data$iso3c, "has no solution for top-to-bottom redistribution."));
                          NaN
                        })
  
  # Stats for the rich population
  wealth.rich <- cubintegrate(xf1, thres.rich, Inf, method = "pcubature", 
                              gini.b=gini.base, avg.b=avg.base, min.b=min.base)   # Manually find x=104.6, which can compensate all sub-DLE pop with the DLE threshold value
  share.rich <- cubintegrate(f1, thres.rich, Inf, method = "pcubature", 
                             gini.b=gini.base, avg.b=avg.base, min.b=min.base)
  mean.rich <- wealth.rich$integral / share.rich$integral
                           
  # Return list
  result = list(share.subDLE = share.subDLE, 
                avg.calculate = avg.calculated, 
                wealth.subDLE = wealth.subDLE, 
                mean.subDLE = mean.subDLE, 
                gap.subDLE = gap.subDLE, 
                mean.gap.subDLE = mean.gap.subDLE,
                thres.rich = thres.rich,
                share.rich = share.rich$integral,
                wealth.rich = wealth.rich$integral,
                mean.rich = mean.rich)
  
  return(list(input = cty.data, result = result))
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


#' @param dle.growth growth constraint for the DLE redistribution case. 
#' @details "max": maximum growth & minimum gini dev
#' @details "no" : no growth 
#' @details "10%" : 10% growth per year & corresponding gini dev 

DeriveIneqStat <- function(cty.data, dle.growth = "max") {
  
  gini.base  = cty.data$input$data$gini.base/100  # WDI has pct values.
  avg.base   = cty.data$input$data$avg.base  
  min.base   = cty.data$input$data$min.base  
  yr.base    = cty.data$input$data$year  
  dle.thres  = cty.data$input$data$dle.thres
  sc         = cty.data$input$sc
  thres.rich = cty.data$result$thres.rich
  
  # Chose a scaler for the DLE counterpart
  if (dle.growth == "no") {
    sc.dle = (avg.base - dle.thres)/(avg.base - min.base) # from Eqn (8)
  } else if (dle.growth == "max") {
    sc.dle = sc[1]    # sc for smallest gini change
  } else { # 10% fixed growth/year
    sc.dle = (1.1^(yr.target-yr.base) * avg.base - dle.thres)/avg.base
  }
  
  # Draw the baseline distribution (to analyze top to bottom redistribution)
  n.draw <- 1e7
  X <- sort(DrawRefLognorm(n.draw, gini.base, avg.base, min.base) + min.base)
  X.dle <- (X - min.base) * sc.dle + dle.thres # Try one example of DLE distributions

  gini.base.draw   = Gini(X)
  gini.dle    = Gini(X.dle)  
  gini.dle.calc = gini.base * (sc.dle * avg.base / (sc.dle * avg.base + dle.thres))  # Eqn (6)
  
  # When there is a solution for thres.rich (i.e. dle.thres is not too high)
  if (!is.nan(thres.rich)) { 
    X.redist <- X
    X.redist[X.redist < dle.thres] <- dle.thres
    X.redist[X.redist > thres.rich] <- thres.rich
    X.redist.jit <- X.redist + runif(n.draw, -1, 1)*1e-4 # Add jitter since lower deciles can have identical values.
    gini.redist = Gini(X.redist.jit)
  }
  
  # Derive extreme decile ratios
  dd.base <- data.frame(X, 
                        grp=cut2(X, quantile(X, seq(0, 1, 0.1))))
  dd.dle <- data.frame(X.dle, 
                       grp=cut2(X.dle, quantile(X.dle, seq(0, 1, 0.1))))
  
  dr.base = dd.base %>% group_by(grp) %>% summarise(s = sum(X)) %>% 
    ungroup() %>% summarise(dr = max(s)/min(s)) %>% as.numeric()
  dr.dle = dd.dle %>% group_by(grp) %>% summarise(s = sum(X.dle)) %>%
    ungroup() %>% summarise(dr = max(s)/min(s)) %>% as.numeric()
  
  if (is.nan(thres.rich)) {
    dr.redist = NaN
    gini.redist = NaN
  } else {
    dd.redist <- data.frame(X.redist.jit, 
                                   grp=cut2(X.redist.jit, quantile(X.redist.jit, seq(0, 1, 0.1)), digits=15))
    dr.redist = dd.redist %>% filter(!is.na(grp)) %>% group_by(grp) %>% summarise(s = sum(X.redist.jit)) %>%
                         ungroup() %>% summarise(dr = max(s)/min(s)) %>% as.numeric()
  }
  
  means = list(mean(X), 
               ifelse(is.nan(thres.rich), NaN, mean(X.redist.jit)), 
               mean(X.dle))
  
  growth.r = (mean(X.dle)/mean(X))^(1/(yr.target-yr.base))
  
  return(list(gini.base      = gini.base, 
              # gini.base.draw = gini.base.draw, 
              gini.redist    = gini.redist, 
              gini.dle       = gini.dle, 
              gini.dle.calc  = gini.dle.calc,
              dr.base        = dr.base, 
              dr.redist      = dr.redist, 
              dr.dle         = dr.dle,
              growth.r       = growth.r,
              means          = means))
  
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

PlotIndiffCurve <- function(cty.data) {
  library(directlabels)
  
  gini.base  = cty.data$input$data$gini.base/100  # WDI has pct values.
  avg.base   = cty.data$input$data$avg.base  
  min.base   = cty.data$input$data$min.base  
  yr.base    = cty.data$input$data$year  
  dle.thres  = cty.data$input$data$dle.thres
  sc         = cty.data$input$sc
  
  # avg.new    = sc * avg.base + dle.thres
  # thres.seq  = seq(1.9*365, dle.thres, 500) # Let's not use 1.9/day threshold
  thres.seq  = seq(dle.thres, 0, -500)[1:4] # dle.thres is based on low-income countries' mean at yr.base.
  
  # growth.r = (avg.new/avg.base)^(1/(yr.target-yr.base))
  # gini.dle.calc = gini.base * (avg.new - dle.thres) / avg.new
  
  l = list(thres = thres.seq, sc = sc)
  df = expand.grid(l) %>% mutate(avg.new = sc * avg.base + thres) %>%
    mutate(growth.r = (avg.new/avg.base)^(1/(yr.target-yr.base)) - 1) %>%
    mutate(gini.dle.calc = gini.base * (avg.new - thres) / avg.new)  
  
  # historical %>% left_join(master.sub) %>% filter(!is.na(gini.base))
  p = ggplot() +
    geom_line(data=df, aes(gini.dle.calc,  growth.r, group = thres)) +
    geom_point(data=historical %>% filter(iso3c == cty.data$input$data$iso3c), 
               aes(gini/100, gr/100, colour = recent), size=3) +
    geom_text(data=historical %>% filter(iso3c == cty.data$input$data$iso3c), 
              aes(gini/100 - 0.01, gr/100, label = year)) +
    labs(title=countrycode(cty.data$input$data$iso3c, 'iso3c', 'country.name')) + scale_x_reverse() +
    labs(x = "Gini index", y = "Annual average growth rate") +
    scale_color_discrete(name = "Year", labels = c("Past", "Latest")) +
    theme(legend.position = c(.95, .95),
          legend.justification = c("right", "top")) +
    # Add threshold value labels
    geom_dl(data=df, aes(gini.dle.calc,  growth.r, 
                         label = paste0('$', format(thres/365, digits=2), '/day')), 
            method = list("last.points", dl.trans(x=x-1.2, y=y+0.5)), cex = 0.7)
  
  return(p)
}


AddRedistLine <- function(p, ineq) {
  gini.base   = ineq$gini.base  # WDI has pct values.
  gini.redist = ineq$gini.redist  # WDI has pct values.
  
  if (!is.nan(gini.redist)) {
    p = p + geom_segment(aes(x = gini.base, y = 0, xend = gini.redist, yend = 0), 
                         arrow = arrow(length = unit(0.03, "npc"))) +
      geom_point(data= data.frame(x = c(gini.base, gini.redist), y=0), aes(x, y), size=3)
  }
  
  return(p)
}

