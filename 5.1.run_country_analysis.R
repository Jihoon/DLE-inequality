
source("5.country_analysis.R")

yr.target <- 2030 #2050 #
yr.base <- 2016

historical <- WDI(country = c("IN", "NE", "RW", "ZA"), indicator = c("NY.GDP.PCAP.KD.ZG", "SI.POV.GINI")) %>% 
  mutate(iso3c = countrycode(iso2c, 'iso2c', 'iso3c'))  %>% rename(gr = NY.GDP.PCAP.KD.ZG, gini = SI.POV.GINI) %>% 
  filter(year <= yr.base) %>%
  group_by(iso3c) 

historical <- historical %>%
  filter(!is.na(gini)) %>% mutate(recent = ifelse(year == max(year), "Latest", "Past")) %>%
  rbind(historical %>% filter(iso3c=="IND", year==yr.base) %>% mutate(gini=37.8, recent="Assumed")) %>%
  arrange(country) %>% mutate(recent = factor(recent, levels=c("Latest", "Past", "Assumed")))
  


####

# Harmonize the base year for gdp.pcap to yr.base
gdp.pcap.base <- raw.gdp.pcap %>% filter(iso3c %in% c('RWA', 'NER', 'IND', 'ZAF', 'LIC'), year == yr.base) %>%
  rename(avg.base=GDP.PCAP) %>% select(-region, -country)

# Threshold for average of 'low income' country groups
gdp.thres <- gdp.pcap.base %>% filter(iso3c=="LIC") %>% select(avg.base) %>% as.numeric()

scaler_1.9 <- WDI(country = c("IN", "RW"), indicator = c("NY.GDP.DEFL.KD.ZG"), start=2012, end = yr.base) %>% 
  rename(r=NY.GDP.DEFL.KD.ZG) %>% mutate(r=r/100+1) %>% group_by(country) %>% summarise(r.tot = prod(r)) %>% 
  ungroup() %>% summarise(r = mean(r.tot))

# Threshold for equivalant of 1.9$/day for yr.base (2016)
lowest.thres <- as.numeric(1.9 * scaler_1.9 * 365)



master.sub.wb <- gdp.pcap.base %>% filter(iso3c != "LIC") %>% 
  left_join(historical %>% filter(recent=="Latest") %>% select(iso3c, gini), by="iso3c") %>%
  mutate(min.base=0, dle.thres = gdp.thres) %>%   # This is for plotting IND illustrative lognorm curves.
  # mutate(min.base=0, dle.thres = lowest.thres) %>%  # This is for plotting indifference curve sets.
  rename(gini.base = gini) %>% select(iso3c, year, everything())

country.list <- split(master.sub.wb, seq(nrow(master.sub.wb)))
names(country.list) <- master.sub.wb$iso3c

sc.list <- lapply(country.list, GetScaler)
y0.list <- lapply(sc.list, GetRefLognorm)
df.list <- lapply(y0.list, GetDF)
result.list <- lapply(df.list, GetBaselinePDF) # This is where distribution type needs to be specified. Before this, it doesn't know about distribution.
ineq.list <- lapply(result.list, DeriveIneqStat, dle.growth="max") # dle.growth= "10%" | "no" | "max"

p.list <- lapply(result.list, PlotIndiffCurve) 
p.list.redist <- mapply(AddRedistLine, p.list, ineq.list, SIMPLIFY = FALSE) 
p.list.redist[[1]]
p.list.redist[[3]]

ExportPDFPlot <- function(name) {
  pdf(file = paste0("plots/Growth-Gini plot ", name, " ",  yr.target,".pdf"), width = 10, height = 6)
  # print(p.list[[name]])       # Not add the no-growth line (for 2050)
  print(p.list.redist[[name]])  # Add the no-growth line (for 2030)
  dev.off()
}
sapply(names(p.list.redist), ExportPDFPlot)
