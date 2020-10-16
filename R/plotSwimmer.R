###########
# plotSwimmer.R
# Created by: Cindy Yang
# Updated: Oct 16, 2020
# This is a function to create Figure 2 swimmer plot in for Bratman, Yang et al. 2020 Nature Cancer
#
# Parameters
# -----------
# plot.dat : dataframe, required 
#      data-frame of plotting data loaded into environment 
#      patientIDs as rowname 
# pt.select : string, optional 
#      string of patientIDs in final plotting order
#      plotting will be in order of plot.dat if NA
# main.lab: string, optional, default = ""
#      string of plot title
# border: boolean, default = FALSE
#      boolean control whether border is drawn around legend
# plot.legend: boolean, default = TRUE
#      boolean control whether legend is drawn
# decision: boolean, default = FALSE (depreciated)
#
# Returns
# ------------
# list: parameters to to regenerate plot
#    list(plot.dat, mp, ctDNAsummary,n.bar, bar.gap, bar.width, y.max, y.lab)
###########

plotSwimmer <- function (pt.select, 
						 main.lab = "", 
						 decision = FALSE,
                         border = FALSE, 
                         plot.legend = TRUE){
  # 0. Set bar colors
  library(RColorBrewer)
  barcol <- c(brewer.pal(12, "Paired")[c(11,5,6,1)], "grey")
  barcol <- c("#FDD901",brewer.pal(12, "Paired")[c(5,6,1)], "grey")
  
  # 1. Subset patients for plotting
  if (!is.na(pt.select)){
  	library(plyr)
    plot.dat <- plot.dat[match(pt.select,row.names(plot.dat)),]
    mp <- barplot (plot.dat[,"Max.time"], plot= FALSE)[,1]
    timeline <- timeline[as.character(timeline$PID) %in% pt.select,]
    timeline$Y.pos <- as.numeric(factor(as.character(timeline$PID), 
                                        levels = row.names(plot.dat)))
    print(head(timeline))
    ctDNA.plot <- ctDNA.plot[as.character(ctDNA.plot$PATIENT_ID) %in% pt.select,]
    ctDNA.plot$Y.pos <- mp[as.numeric(factor(as.character(ctDNA.plot$PATIENT_ID),
                                             levels = row.names(plot.dat)))]
    max.length = ddply (ctDNA.plot, ~PATIENT_ID, summarize, max.val = max(START_DATE_weeks, na.rm = TRUE))
    max.length$Y.pos= mp[as.numeric(factor(as.character(max.length$PATIENT_ID), 
                                           levels = row.names(plot.dat)))]
  }
  
  # 2. Plotting start here
  # Plot axis labels
  x.lab <- "Time since initiation of therapy, wks"
  y.lab <- "Subjects recieved study drug, pts"
  x.lim <- c(0,round_any(as.numeric( max(plot.dat$Max.time, na.rm = TRUE)+10), 20, f = ceiling))
  print(x.lim)
  opar <- par(lwd = 0.3)
  
  # 2.1 Calculate variables that are re-used by multiple functions
  n.bar <- nrow(plot.dat)
  bar.gap <- 0.2
  bar.width <- 1
  y.max <- n.bar*bar.width + (1+n.bar)*bar.gap
  
  # 2.2 Setup empty plot area
  par(mar = c(6.15 ,5, 5.2, 2.1))
  plot(NA, 
       ylim = c(0,y.max), 
       xlim = x.lim, 
       xlab = "", 
       ylab = "",
       yaxs="i",
       xaxt = "n", 
       yaxt = "n")
       
  # 2.3 Plot each line segment
  segments(x0 = timeline$interval.start, 
           x1 = timeline$interval.end, 
           y0 = mp[timeline$Y.pos],
           col = timeline$colors,
           lwd = 5)
 
  # 2.4 Plot location of end of trial
  points (
    x = plot.dat[plot.dat[,"Still on treatment (y or n)"]=="Y","Max.time"]+1.5,
    y = mp[plot.dat[,"Still on treatment (y or n)"]=="Y"],
    pch = -c(8594),
    cex = 1,
    col = "black")
  
  # 2.5 Plot arrow for those continuing treatment
  points (
    x = plot.dat[plot.dat[,"Still on treatment (y or n)"]=="N","Weeks.on.treatment"],
    y = mp[plot.dat[,"Still on treatment (y or n)"]=="N"],
    pch = 23,
    cex = 0.9,
    lwd = 1,
    col = "black")
  
  ## 2.6 Plot symbol for those off trial not due to progression
  #if(decision){
  #  reasons.keep <- "Death"
  #  points (
  #    x = plot.dat[plot.dat[,"REASON_OFF_TRIAL"] %in% reasons.keep,"Weeks.on.treatment"],
  #    y = mp[plot.dat[,"REASON_OFF_TRIAL"] %in% reasons.keep],
  #    pch = 4,
  #    lwd = 2,
  #    cex = 0.70,
  #    col = "black")
  #}
  
  # 2.7 Plot deaths
  points (
    x = plot.dat[as.character(plot.dat[,"OS_STATUS"]) %in% "DECEASED","Max.time"],
    y = mp[as.character(plot.dat[,"OS_STATUS"]) %in% "DECEASED"],
    pch = 4,
    lwd = 2,
    cex = 0.70,
    col = "black") 
                                                       
  # 2.8 ctDNA detection status
  points(
    x = ifelse(ctDNA.plot$START_DATE_weeks<0, 0, ctDNA.plot$START_DATE_weeks),
    y = ctDNA.plot$Y.pos,
    col = "black",
    bg = as.character(ctDNA.plot$DetectedStatus),
    pch = 21,
    cex = 0.70
  )
  
  # 2.9 Add Plot title, axes labels and frame
  title(xlab = x.lab, main = "RECIST 1.1 Clincal Response", cex.lab = 0.7, line = 2.50, cex.main = 0.7)
  axis (side = 1, lwd = 1, labels = T, lwd.tick=1, xaxs="i", cex.axis = 0.7)
  axis (side = 2, lwd = 1, labels = row.names(plot.dat), at = mp, lwd.tick=1, xaxs="i", cex.axis = 0.5, las = 2)
  box(lwd = 0.75)
  
  # 3. Add Figure legend
  swimmer.legend = list(
    labels = c("Stable disease",
               "Partial response",
               "Complete response",
               "Disease progression",
               "No response data",
               "Continued response",
               "Death", 
               "End of trial"),
    ltys = c(1,1,1,1,1, NA, NA, NA),
    lwds = c(2.5,2.5,2.5,2.5,2.5,NA,1.5, 1.2),
    pchs = c(NA,NA,NA,NA,NA,-8594,4, 23),
    pt.cexs = c(NA,NA,NA,NA,NA,1.5,1.3, 1),
    pt.cols = c(barcol,"black","black", "black")
  );
  
  coord <- legend("bottomright",
                  title = "RECIST 1.1",
                  legend = swimmer.legend$labels,
                  pch = swimmer.legend$pchs,
                  lty = swimmer.legend$ltys,
                  col = swimmer.legend$pt.cols,
                  pt.cex = swimmer.legend$pt.cexs,
                  lwd = swimmer.legend$lwds,
                  bty = "n", plot = FALSE,
                  cex = 0.70)
  
  if (plot.legend){
    legend(x= coord$rect$left, y = coord$rect$top+5,
           title = "RECIST 1.1",
           legend = swimmer.legend$labels,
           pch = swimmer.legend$pchs,
           lty = swimmer.legend$ltys,
           col = swimmer.legend$pt.cols,
           pt.cex = swimmer.legend$pt.cexs,
           lwd = swimmer.legend$lwds,
           bty = "n",
           cex = 0.70)
    
    legend(x= coord$rect$left+4, y = coord$rect$top+10,
           title = "ctDNA",
           legend = c("Detected", "Not-detected"),
           pch = c(16,21),
           bg = c("black", "white"),
           pt.cex = c(1.5, 1.5),
           bty = "n",
           cex = 0.70)
  }
  
   # 4. Return plot data
  return (list(plot.dat = plot.dat, mp=mp,
               ctDNAsummary = ctDNAsummary,
               n.bar =n.bar, bar.gap = bar.gap, bar.width = bar.width, y.max = y.max, y.lab = y.lab)
               )
}