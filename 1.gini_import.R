library(readxl)
library(tidyverse)
library(WDI)
library(ggplot2)
library(cubature)
library(Hmisc)
require(reldist)
require(ineq)
require(stats)
require(EnvStats)
library(stats)
library(devtools)
library(countrycode)

source("3.functions.R")

### Data import

# Gini

# 1. Online source for Gini
# raw.gini <- read_xlsx("_Gini Data - v3 - by Gapminder.xlsx", sheet=2) %>% filter(time > 1950) %>% mutate(iso3c=toupper(geo)) %>% select(-geo) %>%
#   rename(year=time)
# 2. From NR (Some values and countries are deleted from orginal data.)
raw.gini <- read.csv("Ginis_wiid34_cleaned.csv") %>% filter(!is.na(rmultiGini)) %>% arrange(Countrycode3, year) %>% 
  rename(iso3c=Countrycode3, Gini.org=Gini, Gini=rmultiGini, name=Country) %>% mutate(diff=Gini-Gini.org) %>% filter(iso3c!="AZE") 
# Aux: WIID 
raw.wiid <- read_xlsx("WIID_19Dec2018.xlsx") %>% filter(source_comments == "PovcalNet") %>%
  rename(iso3c=c3,  Gini=gini_reported, Country=country) %>% 
  select(iso3c, year, Country, Gini, starts_with('source')) %>% arrange(iso3c, -year)

# GDP per cap, PPP (constant 2017 international $)
raw.gdp.pcap <- WDI(country = "all", indicator = "NY.GDP.PCAP.PP.KD", 
                    start = 1990, end = 2019, extra = TRUE, cache = NULL) %>%
  # filter(region!="Aggregates") %>% 
  select(-iso2c, -(capital:lending)) 
names(raw.gdp.pcap)[2] <- 'GDP.PCAP'

raw.pop <- WDI(country = "all", indicator = "SP.POP.TOTL", 
               start = 1990, end = 2019, extra = TRUE, cache = NULL) %>%
  filter(region!="Aggregates") %>% select(-iso2c, -(capital:lending)) 
names(raw.pop)[2] <- 'Population'


### Construct master DF
master <- raw.gdp.pcap %>% left_join(raw.gini) %>% arrange(iso3c, -year) %>% 
  select(country, iso3c, region, year, Gini, GDP.PCAP, -name) %>% 
  filter(!is.na(Gini)) %>% group_by(iso3c) %>%
  mutate(sd=sd(Gini), min.Gini=min(Gini), max.Gini=max(Gini)) %>% 
  mutate(min.year=(Gini==min.Gini), max.year=(Gini==max.Gini)) %>% 
  mutate(ratio=min.Gini/max.Gini) %>% filter(sd>0) 


master %>% filter(min.year | max.year)
master.region <- master %>% group_by(region)  %>% 
  filter(ratio==min(ratio))# & region!="North America")

master.improve <- master %>% filter(first(Gini) < last(Gini))
a <- lm(Gini~year + region + year:region, data=master.improve)
summary(a)

GDP.avg.by.region <- master.improve %>% left_join(raw.pop) %>% group_by(region, year) %>% 
  summarise(avg.gdp.pcap=weighted.mean(GDP.PCAP, Population, na.rm=TRUE))

### Plot - Countries with largest changes in each region
p <- ggplot(master.region, aes(x=year, group=country, label=country))
p + 
  geom_line(aes(y=Gini, color="Gini"), size=1) + 
  geom_line(aes(y=GDP.PCAP/300, color="GDP.PCAP"), size=1, linetype="longdash") +
  geom_label(x=1995, y=60, size=4) +
  facet_wrap(~region) + 
  labs(colour = "Parameter") +
  scale_y_continuous(sec.axis = sec_axis(~.*300, name = "GDP.PCAP [PPP 2011$]"))


### Plot - Scatter dots in each region
p <- ggplot(master, aes(x=jitter(year, amount=0.5), y=Gini, group=region, label=region))
p +
  geom_point(aes(color=region))+
  facet_wrap(~region)  + theme(legend.position = "none") 
  

p <- ggplot(master.improve, aes(x=jitter(year, amount=0.5), y=Gini, group=region, label=region))
p +
  geom_point(aes(color=region))+
  stat_smooth(method = "lm", col = "red") +
  geom_line(data=GDP.avg.by.region, aes(x=year, y=avg.gdp.pcap/700, group=region)) +
  facet_wrap(~region)  + theme(legend.position = "none") 
  

### Getting moving slopes for gini
library(zoo)

slop_e <- function(z) coef(lm(Gini ~ year, as.data.frame(z)))[[2]] # Not used (meant to be for rollapply)
slop_10yr <- function(z) coef(lm(Gini ~ year, z2[z,]))[[2]]

# a <- data.frame(year=1951:2015)
# a <- expand.grid(year=1951:2015, name=unique(raw.gini$name)) 
# a <- raw.gini %>% right_join(a) %>% arrange(name, year)
# 
# z2 <- a %>% select(year, Gini)

z2 <- raw.gini %>% left_join(unique(raw.gdp.pcap %>% select(iso3c, region)), by="iso3c") %>% select(-Gini.org, -r1Gini, -diff)
country <- unique(z2$name)
window.size <- 10 # 10 year moving average

# Derive a column of offset values to mark the end of the 10-year window (for rolling slope), i.e. i:(i+window10[i]) is the window for each obs.
window <- function(cty) { # for each country
  z <- z2 %>% filter(name==cty)
  w <- z$year + window.size
  z$window10 <- sapply(w, function(x) {max(which(z$year<x))})
  # This will give Inf slope for countries with 1< n.obs <window.size. For them let's have just one average slope.
  z <- z %>% mutate_cond(window10=0, dim(z)[1]==window10 & year>max(year)-window.size+1) 
  z <- z %>% mutate_cond(window10=window10-row_number(), window10>0)
  return(z)
  }

# Apply 'window' for all countries to get offsets and accumulate as a df 
Gini.slope.master <- lapply(country, window)
Gini.slope.master <- do.call(rbind, Gini.slope.master)

# A list having all offset values for corresponding rows
w.offset <- lapply(Gini.slope.master$window10, function(x) {seq(0,x)})

# A list having all index values specifying rolling windows for corresponding rows
w <- lapply(1:dim(raw.gini)[1], function(x) {x+w.offset[[x]]})

# 10-yr rolling slope values (%p per 10 years)
Gini.slope.master$slope10 <- sapply(w, slop_10yr)*window.size

# Number of obs used for slope calc
Gini.slope.master <- Gini.slope.master %>% mutate(n.obs = ifelse(sapply(w, length)>1, sapply(w, length), 0)) %>% group_by(name) %>% 
  mutate(min.Gini=min(Gini), min.slope=min(slope10, na.rm = TRUE)) #%>% mutate(n.obs.min.slope=n.obs[slope10==min.slope])

country.inf <- unique(Gini.slope.master$name[is.infinite(Gini.slope.master$min.slope)])

country.inf <- data.frame(name=country.inf, 
                          slope10=window.size*sapply(country.inf, function(cty) coef(lm(Gini ~ year, Gini.slope.master %>% filter(name==cty)))[[2]]))

Gini.slope.master <- Gini.slope.master %>% left_join(country.inf, by="name") %>% 
  mutate(slope10=ifelse(is.na(slope10.x), slope10.y, slope10.x)) %>% #group_by(name) %>%
  add_tally() %>% 
  mutate_cond(n.obs=n, is.infinite(min.slope)) %>%
  mutate_cond(min.slope=slope10, is.infinite(min.slope)) %>% select(-slope10.x, -slope10.y)

Gini.slope.summary <- Gini.slope.master %>% arrange(name, slope10) %>% summarise_all(first) %>% 
  select(iso3c, name, region, min.Gini, min.slope, n.obs, init.Gini=Gini, init.year=year) # Keeping the Gini and year for the period with the maximum slope.

# Show the top-3 countries in each region
Gini.slope.plot <- Gini.slope.summary %>% filter(n.obs>3) %>% 
  # group_by(region) %>% 
  arrange(min.slope) %>% slice(1:15) %>% mutate(period=paste0(init.year, "-", init.year+9))
view(Gini.slope.plot)
# Show the top-3 countries under certain starting Gini
View(Gini.slope.summary %>% filter(init.Gini<=40 & n.obs>3) %>% group_by(region) %>% arrange(region, min.slope) %>% slice(1:3))




lm(slope10 ~ Gini, data=Gini.slope.master)
ggplot(Gini.slope.master %>% filter(n.obs>3 & slope10 <20 & year>=1970), aes(x=jitter(year, amount=0.5), y=slope10, group=region, label=region)) + 
  geom_point(aes(color=region)) + 
  # stat_smooth(method = "lm", col = "red") +
  facet_wrap(~region)  #+ theme(legend.position = "none") 
ggplot(Gini.slope.master %>% filter(n.obs>3 & slope10 <20), aes(x=Gini, y=slope10, group=region, label=region)) + 
  geom_point(aes(color=region)) + 
  # stat_smooth(method = "lm", col = "red") +
  facet_wrap(~region)  #+ theme(legend.position = "none") 

# describe(Gini.slope.master$slope10) # 95% of all obs above -8.2%p/10yr
# describe(Gini.slope.master$Gini) # 95% of all obs above -8.2%p/10yr
quantile((Gini.slope.master %>% filter(n.obs>3))$slope10, probs=seq(0,1,0.01), na.rm=TRUE)

ggplot(raw.gini %>% filter(name=="Slovakia"), aes(x=year, y=Gini)) + geom_point()

# rollapply is a failure with variable window sizes. Also with the full df (z2), slope calculation is wrong.
# z1 <- z2 %>% select(year, Gini)
# z2$slope10_1 <- rollapply(data=z1, width=w.offset, FUN=slop_e, by.column = FALSE, fill = NA, align = "left") *10  # This cannot fine 'Gini' object. (Why?)
# z2$slope10_2 <- rollapply(data=z2, width=2, FUN=slop_e, by.column = FALSE, fill = NA, align = "left") *10  # Doing something wrong with z2 (no idea)




  
