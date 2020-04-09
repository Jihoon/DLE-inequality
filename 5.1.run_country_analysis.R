
source("5.country_analysis.R")

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
ineq.list <- lapply(result.list, DeriveIneqStat, dle.growth="max") # dle.growth= "10%" | "no" | "max"
