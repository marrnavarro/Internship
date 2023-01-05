## Script to generate a volcano plot and a histogram with the p-values from the regression coefficients ##

## Libraries
library(ggplot2)
library(ggrepel)

## Load data
#I used this script for more than one file but I load here only one for simplification
coef<-get(load("Z:/Analysis_projects/MN_analyses/mixed_effects/rbase_scaling/new/fatigue_compound_coef.RDATA"))
coef<-coef$compound

## Volcano plot
ggplot(coef, aes(x = Estimate, y = -log10(p.val)))+
  geom_point(size= 1, col = ifelse(coef$FDR < 0.05, "black","grey")) +
  geom_label_repel(label = ifelse(coef$ENSG %in% coef$ENSG[order(coef$p.val)][1:4], coef$hgnc_symbol, ""), max.overlaps = 100,  force = 25, color = "black", size = 6, segment.size = 0.25, fontface="italic") +
  xlab("Compound Score effect size") +
  ylab("-10log(p-value)") +
  scale_y_continuous(limits = c(0, 7), 
                     breaks = seq(0, 7, 1),
                     sec.axis = dup_axis(breaks = derive(), labels = derive(), name = NULL))+
  scale_x_continuous(limits = c(-3, 3), 
                     breaks = seq(-3, 3, 1)) +
  labs(tag ="A") +
  theme(legend.position = "none", 
        aspect.ratio = 1.5,
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill="white"),
        panel.grid.major.x = element_line(size = 0.25, color = "grey"),
        panel.grid.major.y = element_line(size = 0.25, color = "grey"),
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color ="black", size = 15, face="bold"),
        axis.title.y = element_text(color = "black", size = 15, face="bold"),
        plot.title = element_text(color = "black", size = 18),
        plot.tag = element_text(color ="black", size= 20, face="bold"),
        plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm")
  )

## Histogram
hist(coef$p.val,breaks = 54,main='',xlab = 'P-value')
