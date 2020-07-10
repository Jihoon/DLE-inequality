
source("5.country_analysis.R")

yr.target <- 2030
yr.base <- 2016

historical <- WDI(country = c("IN", "NE", "RW", "ZA"), indicator = c("NY.GDP.PCAP.KD.ZG", "SI.POV.GINI")) %>% 
  mutate(iso3c = countrycode(iso2c, 'iso2c', 'iso3c'))  %>% rename(gr = NY.GDP.PCAP.KD.ZG, gini = SI.POV.GINI) %>% 
  filter(!is.na(gini)) %>% group_by(iso3c) %>% mutate(recent = (year == max(year))) %>% ungroup()

####

# Harmonize the base year for gdp.pcap to yr.base
gdp.pcap.base <- raw.gdp.pcap %>% filter(iso3c %in% c('RWA', 'NER', 'IND', 'ZAF', 'LIC'), year == yr.base) %>%
  rename(avg.base=GDP.PCAP) %>% select(-region, -country)

gdp.thres <- gdp.pcap.base %>% filter(iso3c=="LIC") %>% select(avg.base) %>% as.numeric()

# master.sub <- master %>% filter(iso3c %in% c('RWA', 'NER', 'IND', 'ZAF')) %>%
#   group_by(iso3c) %>% 
#   summarise_all(first) # Keep only the most recent obs.
# 
# master.sub <- master.sub %>% rename(gini.base = Gini, avg.base = GDP.PCAP) %>%
#   # mutate(min.base = 0, dle.thres = 1.9*365) %>% 
#   mutate(min.base = 0, dle.thres = gdp.thres) %>% # mean of 'low-income' countries https://data.worldbank.org/indicator/NY.GDP.PCAP.PP.KD?locations=XM
#   select(-(sd:ratio), -country, -region)

# 02 Jul 2020. We decide to use only WB values for consistency.

# master.sub.wb <- master.sub %>% full_join(historical %>% filter(recent)) %>%
#   mutate(gini.base = coalesce(gini, gini.base)) %>% 
#   left_join(raw.gdp.pcap) %>%
#   mutate(avg.base = coalesce(avg.base, GDP.PCAP), min.base = 0, dle.thres = gdp.thres) %>% # source given above
#   select(names(master.sub)) %>% group_by(iso3c) %>%
#   filter(year == max(year))

master.sub.wb <- gdp.pcap.base %>% filter(iso3c != "LIC") %>% 
  left_join(historical%>% filter(recent) %>% select(iso3c, gini), by="iso3c") %>%
  mutate(min.base=0, dle.thres = gdp.thres) %>%
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
ExportPDFPlot <- function(name) {
  pdf(file = paste0("plots/Growth-Gini plot ", name, ".pdf"), width = 10, height = 6)
  print(p.list.redist[[name]])
  dev.off()
}
sapply(names(p.list.redist), ExportPDFPlot)
