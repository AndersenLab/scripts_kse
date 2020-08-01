library(qtlDesign)
dd = seq(0.01,10,.01)

n100 = power.t.test(n=100, delta=dd, sig.level=.05/750)$power #old WI panel

pvar= prop.var('ri', dd/2,1)

df <- data.frame(pvar, n100)

ggplot(data=df) + aes(x=pvar*100, y=n100*100) + 
  geom_line(size=0.5, color="black", linetype=2) + 
  geom_hline(yintercept = 80) +
  geom_vline(xintercept = 11) +
  xlab("% Variance Exp.") + ylab("Power") + xlim(0, 35) +
  theme(legend.position="none", 
        axis.title.x = element_text(vjust=0, size=12, face="bold"), 
        axis.title.y = element_text(size=12, angle=90, vjust=0.25, face="bold"), 
        axis.text.x=element_text(size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"))
