library(readxl)
library(tidyverse)
library(ggplot2)

input_path = "C:/Users/min/IIASA/DLE - DLE and inequality/data/WIID 2022/"

f = list.files(input_path)
setwd(input_path)

data_raw = readxl::read_xlsx(f[2]) %>% filter(resource != "Earnings")
country_list = c("USA", "BRA", "IND", "IDN", "ZAF", "CHN",
                 "RUS", "RWA")
# country_list = c("IND", "RWA") "DEU", "THA", "CAN", "GBR", 

add_suffix <- function(x) paste0(x, "_share")
change_suffix <- function(x) gsub("_share_pcap", "_pcap", x)

data_filtered = data_raw %>% 
  filter(!is.na(d1), areacovr=="All", #areacovr_detailed=="All", 
         popcovr=="All") %>% #reference_unit=="Household", 
         # population > 1e7) %>%
  filter(c3 %in% country_list) 

data_w = data_filtered %>%
  rename_with(add_suffix, d1:d10) %>%
  mutate(across(d1_share:d10_share, ~ .x*mean/10, .names = "{.col}_pcap")) %>%
  rename_with(change_suffix, d1_share_pcap:d10_share_pcap)

data_l = data_w %>% select(country, c3, gini, region_un, areacovr, areacovr_detailed, reference_unit,
                         year, resource, scale, scale_detailed, source, source_detailed, 
                         survey, d1_share:d10_share, d1_pcap:d10_pcap) %>% 
  pivot_longer(cols=c(d1_share:d10_share, d1_pcap:d10_pcap), names_to = c("decile", "type"), names_sep  = "_") %>%
  pivot_wider(names_from="type") %>%
  mutate(decile = factor(decile, levels = c("d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9", "d10"))) %>%
  group_by(country, c3, year, region_un, areacovr, areacovr_detailed, reference_unit,
           resource, scale, scale_detailed, source, source_detailed, 
           survey) %>%
  mutate(cum_share = cumsum(share))

cty_count = data_filtered %>%
  group_by(country, c3, region_un, areacovr, areacovr_detailed, reference_unit,
           resource, scale, scale_detailed, source, source_detailed) %>% 
           # survey) %>%
  count(country) %>% filter(n>1) %>% arrange(country, -n) %>%
  group_by(country) %>% slice(1)

df_select = cty_count %>% left_join(data_l)

# Lorenz curves by country
ggplot(df_out %>% select(country, year, decile, share, cum_share, cum_share_hat)) +
  geom_line(aes(x=decile, y=share, group=year, colour=year)) +
  # geom_line(aes(x=decile, y=cum_share, group=year, colour=year)) +
  # geom_line(aes(x=decile, y=cum_share_hat, group=year, colour=year)) +
  scale_colour_gradient2(midpoint = 1995, 
                         low = "blue",
                         mid = "yellow",
                         high = "red") +
  # facet_grid(~country, scales="free_x")
  facet_wrap(~country, ncol=4, scales="free_x") +
  # labs(y = 'Cumulative shares of income/consumption [%]', x='Decile')
  labs(y = 'Decile shares among total income/consumption [%]', x='Decile')

# # One country
# ggplot(df_select %>% select(country, year, decile, cum_share, cum_share_hat, pcap) %>%
#              filter(country == "Indonesia", year %in% c(1984, 2015)),
#        aes(x=decile, group=year, colour=year)) +
#   geom_line(aes(y=cum_share, size=year)) +
#   geom_line(aes(y=cum_share_hat), linetype = "dashed") +
#   scale_colour_gradient2(midpoint = 1999.5, mid = "green") 

ggplot(df_select %>% select(country, year, gini)) +
  geom_line(aes(x=year, y=gini)) +
  facet_grid(~country, scales="free_x")

#### quadratic optimization 
# Define x1, x2 as the affine parameters -> c1_hat = c0x1 + x2
# Then find x1 and x2 to minimize OBJ = (c1 - c0x1 - x2)•(c1 - c0x1 - x2) = (c1 - M•x)•(c1 - M•x)
# OBJ = t(c1) %*% c1 - 2*t(c1)•M•x + t(x)•t(M)•M•x, s.t. x1 >= 0, 100x1+x2=100, c01x1+x2 >= 0

library(quadprog)
library(lsa)
library(ineq)

solve_qp = function(c1, c0) {
  
  M = as.matrix(data.frame(a = c0, b=1))
  Dmat = 2 * t(M) %*% M
  dvec = 2 * t(c1) %*% M
  
  # Amat = matrix(c(1,0,0,0), nrow=2)
  # bvec = c(0, 0)
  # qp <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 0)
  
  Amat = matrix(c(100,1, 1,0), nrow=2)
  bvec = c(100, 0)
  # Amat = matrix(c(100,1, c0[1],1, 1,0), nrow=2)
  # bvec = c(100, 0, 0)
  qp <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  
  print(qp$value)
  print(qp$solution)
  return(qp$solution)
}

euclidean <- function(a, b) sqrt(sum((a - b)^2))

df_fy = df_select %>% select(year) %>% slice(which.min(abs(year - median(year)))) %>% rename(first_year = year)
# df_fy = df_select %>% select(year) %>% arrange(year) %>% slice(1) %>% rename(first_year = year)
df_temp = df_select %>% left_join(df_fy) %>% 
  mutate(first_year = (year == first_year)) %>%
  # mutate(first_year = (year == 2005)) %>%
  select(country, c3, gini, resource, year, decile, share, cum_share, first_year)
df_main = df_temp %>% 
  left_join(df_temp %>% filter(first_year) %>% 
              select(decile, share, cum_share) %>% 
              rename(cum_share_base = cum_share, share_base = share))

# df_IDN = a %>% filter(c3 == "IDN") %>% ungroup()

df_param = df_main %>% group_by(country, year) %>%
  summarise(a = solve_qp(cum_share, cum_share_base)[1], b = solve_qp(cum_share, cum_share_base)[2]) %>%
  mutate(share_const = 100*a + b)

df_out = df_main %>% left_join(df_param) %>%
  mutate(cum_share_hat = a * cum_share_base + b) %>% #, share_hat = a * share_base + b) %>%
  group_by(country, year) %>%
  mutate(cum_share_hat_lag = lag(cum_share_hat, default = 0)) %>%
  mutate(share_hat = cum_share_hat - cum_share_hat_lag) %>%
  mutate(similarity = lsa::cosine(cum_share, cum_share_hat), 
         var_r = var(cum_share_hat-cum_share) / var(cum_share),
         dist_r = euclidean(cum_share_hat, cum_share) / euclidean(cum_share, 0)) %>%
  mutate(gini_hat = Gini(share_hat), gini_repl = Gini(share), gini_base = Gini(share_base)) %>%
  mutate(theil_hat = Theil(share_hat), theil_repl = Theil(share), theil_base = Theil(share_base)) %>%
  mutate(gini_affine = gini_base * (mean(share_hat) - b) / mean(share_hat)) %>%
  mutate(palma_hat = share_hat[10]/sum(share_hat[1:4]), palma = share[10]/sum(share[1:4]))
  
df_out %>% summarise(sum_share_hat = sum(share_hat), sum_share = sum(share))

ggplot(df_out %>% select(country, year, similarity), aes(x=year)) +
  geom_line(aes(y=similarity)) +
  facet_grid(~country, scales="free_x")
ggplot(df_out %>% select(country, year, var_r), aes(x=year)) +
  geom_line(aes(y=var_r)) +
  facet_grid(~country, scales="free_x")
ggplot(df_out %>% select(country, year, dist_r), aes(x=year)) +
  geom_line(aes(y=dist_r)) +
  facet_grid(~country, scales="free_x")


df_plot = df_out %>% 
  select(country, year, gini, gini_hat, gini_repl, dist_r, 
         palma_hat, palma, theil_hat, theil_repl) %>%
  left_join(cty_count %>% select(country, n)) %>%
  mutate(cty.lab = paste0(country, " (", n, ")"))
country.labs = unique(df_plot$cty.lab)
names(country.labs) <- unique(df_plot$country)

coeff = 8
ggplot(df_plot, 
       aes(x=year)) +
  # geom_line(aes(y=gini/100), color="green") +
  geom_line(aes(y=gini_hat, color="Estimated", linetype="Gini")) +
  geom_line(aes(y=gini_repl, color="Observed", linetype="Gini")) +
  # geom_line(aes(y=theil_hat/coeff), color="blue", linetype="dashed") +
  # geom_line(aes(y=theil_repl/coeff), color="red", linetype="dashed") +
  geom_line(aes(y=palma_hat/coeff, color="Estimated", linetype="Palma")) +
  geom_line(aes(y=palma/coeff, color="Observed", linetype="Palma")) +
  # geom_line(aes(y=dist_r*coeff), color="darkgreen") +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Gini",
    breaks=seq(0,1,0.2),
    
    # Add a second axis and specify its features
    # sec.axis = sec_axis(~./coeff, name="Euclidean dist")
    # sec.axis = sec_axis(~.*coeff, name="Theil")
    sec.axis = sec_axis(~.*coeff, name="Palma")
  ) +
  facet_wrap(~country, ncol=4, scales="free_x", labeller = labeller(country = country.labs)) +
  scale_color_manual(name = "Type", values = c("Observed" = "red", "Estimated" = "blue")) +
  scale_linetype_manual(name = "Metric", values = c("Gini" = "solid", "Palma" = "dashed"))
  # facet_grid(~country, scales="free_x")






# # Testing/drafting quadprog
# 
# # c0 = df_IDN %>% filter(year == min(df_IDN$year)) %>% select(share) %>% as.matrix()
# # c1 = df_IDN %>% filter(year == 2000) %>% select(share) %>% as.matrix()
# 
# M = as.matrix(data.frame(a = c0, b=1))
# Dmat = 2 * t(M) %*% M
# dvec = 2 * t(c1) %*% M
# # Dmat = 2 * matrix(c(t(c0) %*% c0, 100, 100, 10), nrow = 2)
# # dvec = 2 * matrix(c(t(c1) %*% c0, 10), nrow = 1)
# 
# Amat = matrix(c(1,0,0,0), nrow=2)
# bvec = c(0, 0)
# 
# qp <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 0)
# qp
# val = qp$value + t(c1) %*% c1
# c1_hat = qp$solution[1] * c0 + qp$solution[2]
# c1
# var(c1_hat-c1) / var(c1)
# t(c1)%*%c1_hat / sqrt(t(c1)%*%c1) / sqrt(t(c0)%*%c0) # Cosine similarity
