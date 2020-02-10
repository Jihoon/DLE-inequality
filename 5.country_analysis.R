###
# Countries to analyze for EL: Rwanda, Niger, India, S.Africa = RWA, NER, IND, ZAF
###

master.sub <- master %>% filter(iso3c %in% c('RWA', 'NER', 'IND', 'ZAF')) %>%
  group_by(iso3c) %>% summarise_all(first)
