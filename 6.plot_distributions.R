gtext = data.frame(xv = apply(df[-1], 2, function(x){x0[which(x==max(x, na.rm=T))]}),
                   yv = apply(df[-1], 2, function(x) max(x, na.rm=T)+0.001),
                   c = paste0('x', 0:5),
                   lab = c(paste0("G[x]:", gini.base),
                           paste0("list(G[z]:0.45, r:5.8*'%'/yr)"),
                           paste0("list(G[z]:0.48, r:6.8*'%'/yr)"),
                           paste0("list(G[z]:0.50, r:7.8*'%'/yr)"),
                           paste0("list(G[z]:0.54, r:9.8*'%'/yr)"),
                           paste0("list(G[z]:0.33, r:2.4*'%'/yr)")))


pdf(file = paste0("India GINI illustration-", distr, ".pdf"), width = 10, height = 6)
# Plot for Princeton w/s, May 2019
ggplot(df %>% select(-y2) %>% gather(key=dist, value=f, -x0), aes(x=x0, y=f, color=dist)) +
  geom_line(size=1) +
  geom_vline(xintercept = min.base, linetype="dashed", colour="blue") +
  geom_vline(xintercept = dle.thres, linetype="dashed", colour="brown") +
  geom_text(aes(x=min.base, label=paste0("Actual minimum : ", min.base, "\n"), y=0.02), colour="blue", hjust=0.0, angle=90, size=4.5) +
  geom_text(aes(x=dle.thres, label=paste0("Desired minimum : ", dle.thres, "\n"), y=0.02), colour="brown", hjust=0.0, angle=90, size=4.5) +
  labs(x="Final energy consumption (GJ/cap)", y="Density")+
  # geom_text(aes(x=x0[which(x0==max(x0, na.rm=T))], label=paste0("g_x: ", gini.base),   y=max(x0, na.rm=T)+0.001, colour="x0"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x1==max(x1, na.rm=T))], label=paste0("g_z: 0.45, r: 5.8%"),  y=max(x1, na.rm=T)+0.001, colour="x1"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # # geom_text(aes(x=x0[which(x2==max(x2, na.rm=T))], label=paste0("g: 0.48, r: 6.8%"),  y=max(x2, na.rm=T)+0.001, colour="x2"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x3==max(x3, na.rm=T))], label=paste0("g_z: 0.50, r: 7.8%"),  y=max(x3, na.rm=T)+0.001, colour="x3"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x4==max(x4, na.rm=T))], label=paste0("g_z: 0.54, r: 9.8%"), y=max(x4, na.rm=T)+0.001, colour="x4"), size=4, hjust=0.0, vjust=0.0, data=df) +
  # geom_text(aes(x=x0[which(x5==max(x5, na.rm=T))], label=paste0("g_z: 0.33, r: 2.4%")), y=max(x5, na.rm=T)+0.001, colour="x5"), size=4, hjust=0.0, vjust=0.0, data=df) +
  geom_text(data=gtext[-3,], aes(x=xv, label=lab, y=yv, colour=c), size=4, hjust=0.0, vjust=0.0, parse=TRUE) +
  theme(legend.position = "none") +
  xlim(c(0, 30))

dev.off()
