###########
# plotFigure2.R
# Created by: Cindy Yang
# Updated: Oct 16, 2020
# This script allows the user to to create Figure 2 
# This script requires the plotSwimmer.R function and Data.Figure2.RData 

###########
#####
# SETUP
#####
source("./R/plotSwimmer.R")
load ("./data/Data.Figure2.RData")

#######
# PLOT MAIN
########
# 1. Plot swimmer plot on right hand side
pdf(height = 12, width = 10, file = "./output/Figure2.pdf")
split.screen (
  matrix (c(c(0.12,1,0,1),
            c(0,0.12,0,1)),
    byrow = T, ncol = 4)
);

# 1) Swimmer plot in order based on % change in ctDNA
screen(1);
# get patients with change
pts.selected <- pts.selected[order(pts.selected$ctDNAChangeLog, decreasing = T),];
pts.selected <- pts.selected[!is.na(pts.selected$ctDNAChangeLog),]; # move the patient with NA to the top of the plot order
plot.dat.netera3 <- plotSwimmer(pt.select = c(row.names(pts.selected),"INS-D-012"),decision = TRUE,plot.legend = TRUE);

# 2) left hand side % change in ctDNA 
screen(2);
par(mar = c(6.15 ,1, 5.2, 0));
pts.selected$ctDNAChangeLog[pts.selected$ctDNAChangeLog==Inf] <- 45; # capping maximum at 45
log10change <- log10(abs(pts.selected[,"ctDNAChangeLog"])*100);
log10change <- ifelse (pts.selected[,"ctDNAChangeLog"]<0, yes = -1*log10change, no = log10change);

plot(
  x = c(-1*log10change,NA),
  y = plot.dat.netera3$mp,
  col = adjustcolor(factor(pts.selected$ctDNAChangeLog>0,levels = c(TRUE, FALSE), labels = c("blue", "red")), alpha.f = 1),
  pch = 16,
  cex = 1,
  xlim = c(-5,2.5),
  ylim = c(0,max(plot.dat.netera3$mp)+0.7),
  axes = FALSE,
  xaxs="i",
  yaxs="i",
  xlab = ""
);

axis(1, at=c(-4:2), labels=c(10000, 1000, 100, 10 ,0, -10,-100),lwd = 1, lwd.tick=1, xaxs="i", cex.axis = 0.5);
# Add title and axes
title(xlab = "Percent change mean\nctDNA detected, %", main = "ctDNA Change\nat cycle 3", cex.lab = 0.60, line = 2.50, cex.main = 0.7);

# Add in ND notation for two samples with no ctDNA detected at baseline
text(x = 0, y = max(plot.dat.netera3$mp), labels = "ND", cex = 0.70);
legend("bottomright",
       title = "Change in %\nctDNA mean",
       legend = c("Increase", "Decrease"),
       fill = c("blue", "red"),
       bty = "n",
       cex = 0.60);

# 3. Close all screens
close.screen(all = TRUE);
dev.off()
