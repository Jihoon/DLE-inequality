# gtext = data.frame(xv = apply(df[-1], 2, function(x){x0[which(x==max(x, na.rm=T))]}),
#                    yv = apply(df[-1], 2, function(x) max(x, na.rm=T)+0.001),
#                    c = paste0('x', 0:5),
#                    lab = c(paste0("G[x]:", gini.base),
#                            paste0("list(G[z]:0.45, r:5.8*'%'/yr)"),
#                            paste0("list(G[z]:0.48, r:6.8*'%'/yr)"),
#                            paste0("list(G[z]:0.50, r:7.8*'%'/yr)"),
#                            paste0("list(G[z]:0.54, r:9.8*'%'/yr)"),
#                            paste0("list(G[z]:0.33, r:2.4*'%'/yr)")))


# Data from the analysis from 5.1
df.plot <- df.list$IND$df %>% select(x0, y0, y20, y17, y14, y11, y8) 
df.data <- df.list$IND$data
df.sc <- df.list$IND$sc[c(20, 17, 14, 11, 8)]
names(df.sc) # corresponding gini values
# g = df.data$gini.base * df.sc * df.data$avg.base / (df.sc * df.data$avg.base + df.data$dle.thres)
growth.r <- (((df.sc * df.data$avg.base + df.data$dle.thres) / df.data$avg.base) ^(1/10) - 1)*100
  
gtext = data.frame(xv = apply(df.plot[-1], 2, function(x){df.plot$x0[which(x==max(x, na.rm=T))]}),
                   yv = apply(df.plot[-1], 2, function(x) max(x, na.rm=T)+0.00001),
                   dist = c("y0", "y20", "y17", "y14", "y11", "y8"),
                   lab = c(paste0("G[x]:", format(df.data$gini.base/100, digits=2)),
                           paste0("list(G[z]:", format(as.numeric(names(df.sc)), digits=2), 
                                  ", r:", format(growth.r, digits =1), "*'%'/yr)")))


pdf(file = paste0("India GINI illustration-", distr, ".pdf"), width = 10, height = 6)
# Plot for Princeton w/s, May 2019
df.plot.l <- df.plot %>%  pivot_longer(-x0, names_to="dist") %>%
  mutate(type = ifelse(dist=="y0", "base", "new")) # For line type setting
ggplot(df.plot.l, aes(x=x0, y=value, color=dist)) +
  geom_line(aes(linetype=type), size=1) +
  xlim(c(0, 2e4)) +
  scale_linetype_manual(values=c("longdash", "solid")) +
  geom_vline(xintercept = df.data$min.base, linetype="dashed", colour="blue") +
  geom_vline(xintercept = df.data$dle.thres, linetype="dashed", colour="brown") +
  # geom_text(aes(x=df.data$min.base, label=paste0("Actual minimum : ", df.data$min.base, "\n"), y=5e-4), colour="blue", hjust=0.0, angle=90, size=4.5) +
  geom_text(aes(x=df.data$dle.thres, label=paste0("Desired minimum : $", format(df.data$dle.thres/365, digits=2), "/day", "\n"), y=1e-4), 
            colour="brown", hjust=0.0, angle=90, size=4.5) +
  # labs(x="Final energy consumption (GJ/cap)", y="Density")+
  labs(x="GDP per capita ($/cap)", y="Density")+
  # geom_text(aes(x=x0[which(x0==max(x0, na.rm=T))], label=paste0("g_x: ", gini.base),   y=max(x0, na.rm=T)+0.001, colour="x0"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x1==max(x1, na.rm=T))], label=paste0("g_z: 0.45, r: 5.8%"),  y=max(x1, na.rm=T)+0.001, colour="x1"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # # geom_text(aes(x=x0[which(x2==max(x2, na.rm=T))], label=paste0("g: 0.48, r: 6.8%"),  y=max(x2, na.rm=T)+0.001, colour="x2"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x3==max(x3, na.rm=T))], label=paste0("g_z: 0.50, r: 7.8%"),  y=max(x3, na.rm=T)+0.001, colour="x3"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x4==max(x4, na.rm=T))], label=paste0("g_z: 0.54, r: 9.8%"), y=max(x4, na.rm=T)+0.001, colour="x4"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x5==max(x5, na.rm=T))], label=paste0("g_z: 0.33, r: 2.4%")), y=max(x5, na.rm=T)+0.001, colour="x5"), size=4, hjust=0.0, vjust=0.0, data=df) +
  geom_text(data=gtext, aes(x=xv, label=lab, y=yv, colour=dist), size=4, hjust=0.0, vjust=0.0, parse=TRUE) +
  theme(legend.position = "none") 

dev.off()
