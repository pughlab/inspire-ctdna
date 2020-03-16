###########
# plotFigure2.R
# Created by: Cindy Yang
# Main R script to create Figure 2 
###########
#####
# SETUP
#####
load("~./data/ASCO_SwimmerPlot_ctDNA_combined_April242019.RData")
source("~./R/plotSwimmer.R")
survival.df <- "~./data/data_allpts_updated_20191026_for_cindy.csv"

# Process survival data
survival.df <- read.csv(survival.df, header = T, stringsAsFactors = F)
pts.selected <- data.frame(
  PATIENT_ID = survival.df$INS_No,
  ctDNAChangeLog = survival.df$change,
  row.names = survival.df$INS_No
)

#######
# PLOT MAIN
########
# 1. Plot swimmer plot on right hand side
# get patients with change

pdf(height = 12, width = 10, 
    file = "~./output/Manuscript_Fig_Swimmer_with_change_Oct.pdf")

split.screen (
  matrix (
    c(c(0.250,1,0,1),
      c(0,0.250,0,1)
    ),
    byrow = T, ncol = 4
  )
);

# 1) Swimmer plot in order based on % change in ctDNA
screen(1);
pts.selected <- pts.selected[order(pts.selected$ctDNAChangeLog, decreasing = T),]
pts.selected <- pts.selected[!is.na(pts.selected$ctDNAChangeLog),]
plot.dat.netera3 <- plotSwimmer(pt.select = c(row.names(pts.selected),"INS-D-012"),decision = TRUE, reprocess.data = FALSE, plot.legend = FALSE)

# 2) left hand side % change in ctDNA 
screen(2);
#par(mar =c(7.2,1,4.1,0.1));
par(mar = c(6.15 ,1, 5.2, 0))
pts.selected$ctDNAChangeLog[pts.selected$ctDNAChangeLog==Inf] <- 45
log10change <- log10(abs(pts.selected[match(ordered.pts,pts.selected$PATIENT_ID),"ctDNAChangeLog"])*100)
log10change <- ifelse (pts.selected[match(ordered.pts,pts.selected$PATIENT_ID),"ctDNAChangeLog"]<0,
                       yes = -1*log10change, 
                       no = log10change)
log10change <- log10(abs(pts.selected[,"ctDNAChangeLog"])*100)
log10change <- ifelse (pts.selected[,"ctDNAChangeLog"]<0,
                       yes = -1*log10change, 
                       no = log10change)

plot(
  x = c(-1*log10change,NA),
  y = plot.dat.netera3$mp,
  col = adjustcolor(factor(pts.selected$ctDNAChangeLog>0,levels = c(TRUE, FALSE), labels = c("red", "blue")), alpha.f = 1),
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
title(xlab = "Percent change mean ctDNA detected, %", main = "ctDNA Change at cycle 3", cex.lab = 0.7, line = 2.50, cex.main = 0.7);

# Add in ND notation for two samples with no ctDNA detected at baseline
text(x = 0, y = max(plot.dat.netera3$mp), labels = "ND", cex = 0.70);
legend("bottomright",
       title = "Change in %\nctDNA mean",
       legend = c("Increase", "Decrease"),
       fill = c("red", "blue"),
       bty = "n",
       cex = 0.70);

# 3. Close all screens
close.screen(all = TRUE);
dev.off()
