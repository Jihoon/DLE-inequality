###
# Identify countries to show in the EL script
# Note: This can be run standalone.
###

library(WDI)
library(tidyverse)
# library(countrycode)

cty.pov <- WDI(indicator = c("SI.POV.DDAY", "SI.POV.LMIC", "SI.POV.UMIC", 
                             "NY.GDP.PCAP.PP.CD", "SP.POP.TOTL", "NY.GDP.PCAP.KD.ZG", "NY.GDP.MKTP.KD.ZG"), 
               start = 2010, end = 2017, extra = TRUE, cache = NULL) %>% 
  filter(capital!="") %>%         # Get only countries, no aggregate regions
  filter(!is.na(SI.POV.DDAY))     # only those with nonzero headcounts under poverty
  
cty.pov.cnt <- cty.pov %>% 
  mutate(pop1.9 = SI.POV.DDAY*SP.POP.TOTL,
         pop3.2 = SI.POV.LMIC*SP.POP.TOTL,
         pop5.5 = SI.POV.UMIC*SP.POP.TOTL) %>%
  filter(year==2014) %>%
  left_join(master, by=c("iso3c", "year")) %>% # master has too few countries with GINI information, but a bit more for 2014 than for 2015. (41 vs 30)
  arrange(pop1.9)

cty.pov.selec <- cty.pov %>%
  filter(iso3c %in% c("IND", "ZAF", "MOZ", "NER", "BWA", "RWA")) %>%
  left_join(master, by=c("iso3c", "year")) %>% 
  mutate(pop1.9 = SI.POV.DDAY*SP.POP.TOTL,
         pop3.2 = SI.POV.LMIC*SP.POP.TOTL,
         pop5.5 = SI.POV.UMIC*SP.POP.TOTL)
