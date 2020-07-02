
source("5.country_analysis.R")

yr.target <- 2030

historical <- WDI(country = c("IN", "NE", "RW", "ZA"), indicator = c("NY.GDP.PCAP.KD.ZG", "SI.POV.GINI")) %>% 
  mutate(iso3c = countrycode(iso2c, 'iso2c', 'iso3c'))  %>% rename(gr = NY.GDP.PCAP.KD.ZG, gini = SI.POV.GINI) %>% 
  filter(!is.na(gini)) %>% group_by(iso3c) %>% mutate(recent = (year == max(year))) %>% ungroup()

####
master.sub <- master %>% filter(iso3c %in% c('RWA', 'NER', 'IND', 'ZAF')) %>%
  group_by(iso3c) %>% 
  summarise_all(first) # Keep only the most recent obs.

master.sub <- master.sub %>% rename(gini.base = Gini, avg.base = GDP.PCAP) %>%
  # mutate(min.base = 0, dle.thres = 1.9*365) %>% 
  mutate(min.base = 0, dle.thres = 2287.8) %>% # mean of 'low-income' countries (2018) https://data.worldbank.org/indicator/NY.GDP.PCAP.PP.CD?locations=XM
  select(-(sd:ratio), -country, -region)

country.list <- split(master.sub, seq(nrow(master.sub)))
names(country.list) <- master.sub$iso3c

sc.list <- lapply(country.list, GetScaler)
y0.list <- lapply(sc.list, GetRefLognorm)
df.list <- lapply(y0.list, GetDF)
result.list <- lapply(df.list, GetBaselinePDF) # This is where distribution type needs to be specified. Before this, it doesn't know about distribution.
ineq.list <- lapply(result.list, DeriveIneqStat, dle.growth="max") # dle.growth= "10%" | "no" | "max"

p.list <- lapply(result.list, PlotIndiffCurve) 
p.list.redist <- mapply(AddRedistLine, p.list, ineq.list, SIMPLIFY = FALSE) 
