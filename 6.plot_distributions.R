# Data for India from the analysis from 5.1
# Plot this for more illustrative threshold = gdp.thres
df.plot <- df.list$IND$df %>% select(x0, y0, y107, y101, y96, y92) 
df.data <- df.list$IND$data
df.sc <- df.list$IND$sc[c(107, 101, 96, 92)]
names(df.sc) # corresponding gini values
# g = df.data$gini.base * df.sc * df.data$avg.base / (df.sc * df.data$avg.base + df.data$dle.thres)
growth.r <- (((df.sc * df.data$avg.base + df.data$dle.thres) / df.data$avg.base) ^(1/(yr.target-yr.base)) - 1)*100
  
gtext = data.frame(xv = apply(df.plot[-1], 2, function(x){df.plot$x0[which(x==max(x, na.rm=T))]}),
                   yv = apply(df.plot[-1], 2, function(x) max(x, na.rm=T)+0.00001),
                   dist = c("y0", "y107", "y101", "y96", "y92"),
                   lab = c(paste0("G[x]:", format(df.data$gini.base/100, digits=2)),
                           paste0("list(G[z]:", format(as.numeric(names(df.sc)), digits=2), 
                                  ", r:", format(growth.r, digits =1), "*'%'/yr)")))


# Export to PDF
pdf(file = paste0("plots/India GINI illustration-", distr, " ", yr.target, "b.pdf"), width = 10, height = 6)

df.plot.l <- df.plot %>%  pivot_longer(-x0, names_to="dist") %>%
  mutate(type = ifelse(dist=="y0", "base", "new")) # For line type setting
ggplot(df.plot.l, aes(x=x0, y=value, color=dist)) +
  geom_line(aes(linetype=type), size=1) +
  xlim(c(0, 5e4)) +
  scale_linetype_manual(values=c("longdash", "solid")) +
  geom_vline(xintercept = df.data$min.base, linetype="dashed", colour="blue") +
  geom_vline(xintercept = df.data$dle.thres, linetype="dashed", colour="brown") +
  geom_text(aes(x=df.data$dle.thres, label=paste0("Desired minimum : $", format(df.data$dle.thres/365, digits=2), "/day", "\n"), y=1e-4), 
            colour="brown", hjust=0.0, angle=90, size=4.5) +
  labs(x="GDP per capita ($/cap)", y="Density (population share)")+
  geom_text(data=gtext, aes(x=xv, label=lab, y=yv, colour=dist), size=4, hjust=0.0, vjust=0.0, parse=TRUE) +
  theme(legend.position = "none") 

dev.off()
