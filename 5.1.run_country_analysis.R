
source("5.country_analysis.R")

yr.target <- 2030 #2050 #
yr.base <- 2016

historical <- WDI(country = c("IN", "NE", "RW", "ZA"), indicator = c("NY.GDP.PCAP.KD.ZG", "SI.POV.GINI")) %>% 
  mutate(iso3c = countrycode(iso2c, 'iso2c', 'iso3c'))  %>% rename(gr = NY.GDP.PCAP.KD.ZG, gini = SI.POV.GINI) %>% 
  filter(year <= yr.base) %>%
  group_by(iso3c) 

historical <- historical %>%
  filter(!is.na(gini)) %>%
  rbind(historical %>% filter(iso3c=="IND", year==yr.base) %>% mutate(gini=37.8)) %>%
  mutate(recent = ifelse(year == max(year), "Latest", "Past")) %>%
  arrange(country) # %>% mutate(recent = factor(recent, levels=c("Latest", "Past", "Assumed")))
  


####

# Harmonize the base year for gdp.pcap to yr.base
gdp.pcap.base <- raw.gdp.pcap %>% filter(iso3c %in% c('RWA', 'NER', 'IND', 'ZAF', 'LIC'), year == yr.base) %>%
  rename(avg.base=GDP.PCAP) %>% select(-region, -country)

# Threshold for average of 'low income' country groups
gdp.thres <- gdp.pcap.base %>% filter(iso3c=="LIC") %>% select(avg.base) %>% as.numeric()

scaler_infl <- WDI(country = c("IN", "RW"), indicator = c("NY.GDP.DEFL.KD.ZG"), start=2012, end = yr.base) %>% 
  rename(r=NY.GDP.DEFL.KD.ZG) %>% mutate(r=r/100+1) %>% group_by(country) %>% summarise(r.tot = prod(r)) %>% 
  ungroup() %>% summarise(r = mean(r.tot))

# Household & NPISH share of GDP
sh.hh.NPISH <- WDI(country = c("IN", "NE", "RW", "ZA"), indicator = "NE.CON.PRVT.ZS", start=yr.base, end = yr.base) %>% 
  mutate(iso3c = countrycode(iso2c, 'iso2c', 'iso3c')) %>% mutate(sh = NE.CON.PRVT.ZS/100)

# Threshold for equivalant of 1.9$/day for yr.base (2016)
thres1.val <- as.numeric(1.9 * scaler_infl * 365)
thres2.val <- as.numeric(3.2 * scaler_infl * 365)

# $1.9 is for per-capita income, need to translate it to per-cap GDP
thres1 <- as.list(thres1.val / sh.hh.NPISH$sh)
thres2 <- as.list(thres2.val / sh.hh.NPISH$sh)
names(thres1) <- sh.hh.NPISH$iso3c
names(thres2) <- sh.hh.NPISH$iso3c


master.sub.wb <- gdp.pcap.base %>% filter(iso3c != "LIC") %>% 
  left_join(historical %>% filter(recent=="Latest") %>% select(iso3c, gini), by="iso3c") %>%
  # mutate(min.base=0, dle.thres = as.numeric(thres2)) %>%   # This is for plotting IND illustrative lognorm curves.
  mutate(min.base=0, dle.thres = as.numeric(thres1)) %>%  # This is for plotting indifference curve sets.
  rename(gini.base = gini) %>% select(iso3c, year, everything())

country.list <- split(master.sub.wb, seq(nrow(master.sub.wb)))
names(country.list) <- master.sub.wb$iso3c

sc.list <- lapply(country.list, GetScaler)
y0.list <- lapply(sc.list, GetRefLognorm)
df.list <- lapply(y0.list, GetDF)
result.list <- lapply(df.list, GetBaselinePDF) # This is where distribution type needs to be specified. Before this, it doesn't know about distribution.

# Param choices: dle.growth= "grow" | "no" | "max", dr.type = "Palma" | "D9_D1" | "10to1"
ineq.list <- lapply(result.list, DeriveIneqStat, dle.growth="max", dr.type = "Palma") # Maximum possible growth for the affine case

# Set specific growth rate given by growth.r (annual per-cap growth)
# 5.4% and 4.8% taken from the manuscript (2006-2016 average)
last.gr <- historical %>% group_by(country) %>% summarise(gr = last(gr)) # Lastest GDP p.c growth rate
ineq.IND <- DeriveIneqStat(result.list$IND, dle.growth="grow", dr.type = "Palma", growth.r = last.gr$gr[last.gr$country=="India"]/100)
ineq.IND.no <- DeriveIneqStat(result.list$IND, dle.growth="no", dr.type = "Palma")
ineq.IND.5_4 <- DeriveIneqStat(result.list$IND, dle.growth="grow", dr.type = "Palma", growth.r = 0.054)
ineq.IND.10 <- DeriveIneqStat(result.list$IND, dle.growth="grow", dr.type = "Palma", growth.r = 0.1)
ineq.IND.max <- DeriveIneqStat(result.list$IND, dle.growth="max", dr.type = "Palma")
ineq.RWA <- DeriveIneqStat(result.list$RWA, dle.growth="grow", dr.type = "Palma", growth.r = last.gr$gr[last.gr$country=="Rwanda"]/100)
ineq.RWA.no <- DeriveIneqStat(result.list$RWA, dle.growth="no", dr.type = "Palma")
ineq.RWA.4_8 <- DeriveIneqStat(result.list$RWA, dle.growth="grow", dr.type = "Palma", growth.r = 0.048)
ineq.RWA.10 <- DeriveIneqStat(result.list$RWA, dle.growth="grow", dr.type = "Palma", growth.r = 0.1)
ineq.RWA.max <- DeriveIneqStat(result.list$RWA, dle.growth="max", dr.type = "Palma")

p.list <- lapply(result.list, PlotIndiffCurve) 
p.list.redist <- mapply(AddRedistLine, p.list, ineq.list, SIMPLIFY = FALSE) 
p.list.redist[[1]]
p.list.redist[[3]]

ExportPDFPlot <- function(name) {
  pdf(file = paste0("plots/Growth-Gini plot ", name, " ",  yr.target,"f.pdf"), width = 10, height = 6)
  print(p.list[[name]])       # Not add the no-growth line (for 2050)
  # print(p.list.redist[[name]])  # Add the no-growth line (for 2030)
  dev.off()
}
sapply(names(p.list.redist), ExportPDFPlot)
sapply(names(p.list), ExportPDFPlot)
