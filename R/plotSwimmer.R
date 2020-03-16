###########
# plotSwimmer.R
# Created by: Cindy Yang
# Function to create swimmer plot 
###########

plotSwimmer <- function (pt.select=NA ,main.lab = "", decision = FALSE, covariates = list(),
                         border = FALSE, barcol = FALSE, stacked = TRUE, reprocess.data = FALSE, plot.legend = TRUE){
  library(plyr)
  # if stacked == TRUE, then color bars blue for time on treatment, orange for progression free time off trial
  if (stacked){
    library(RColorBrewer)
    barcol <- c(brewer.pal(12, "Paired")[c(11,5,6,1)], "grey")
    barcol <- c("#FDD901",brewer.pal(12, "Paired")[c(5,6,1)], "grey")
    
  }
  
  # prepare bar information
  if(reprocess.data) {
    tm.all.tps.df$time.interval <- unlist(lapply(unique(tm.all.tps.df$PID), 
                                                 function (x) { sub.dat = tm.all.tps.df[tm.all.tps.df$PID == x,]
                                                 return (difftime(time1 = as.Date(sub.dat$Date),
                                                                  time2 = as.Date(sub.dat[1,"Date"]), 
                                                                  units = "weeks"))}))
    # get overall response string
    # merge consecutive scans with the same response
    timeline <- do.call ("rbind", lapply(unique(as.character(tm.all.tps.df$PID)), 
                                         function (x) { sub.dat = tm.all.tps.df[tm.all.tps.df$PID == x,]
                                         vals = rle(as.vector(sub.dat[,"Overall.Response"]))
                                         interval.start = cumsum(c(1,vals$length[-length(vals$length)]))
                                         interval.end = interval.start+(vals$length)
                                         if (max(interval.end)>nrow(sub.dat)) { interval.end[which.max(interval.end)] = nrow(sub.dat)}
                                         if (length(vals$values)>1){
                                           interval.length = 
                                             c(sub.dat[interval.end, "time.interval"][-1], plot.dat[x,"Max.time"]) - 
                                             sub.dat[interval.start, "time.interval"]
                                         }else {
                                           interval.length = plot.dat[x,"Max.time"]
                                         }
                                         return(
                                           data.frame(
                                             PID = rep (x, length(vals$values)),
                                             OR = vals$values,
                                             interval.length = interval.length, 
                                             interval.start = sub.dat[interval.start, "time.interval"],
                                             interval.end = sub.dat[interval.start,"time.interval"] + interval.length,
                                             Max.time = rep(plot.dat[x,"Max.time"], length(vals$values))
                                           ))}))
    
    # Y position according to increasing order so the highest OS pt will be on the top of the plot
    
    plot.dat <- plot.dat[unique(as.character(timeline$PID)),]
    plot.dat <- plot.dat[order(plot.dat[,"Max.time"], decreasing = FALSE),]
    timeline$colors <- as.character(factor(as.character(timeline$OR), 
                                           levels = c("SD", "PR", "CR", "PD", "NE"), labels = barcol))
    
    # process ctDNA data for plotting
    plot.dat$WEEKS_SINCE_DIAGNOSIS = difftime(plot.dat$DATE_OF_FIRST_DOSE,
                                              plot.dat$DATE_OF_DIAGNOSIS, units = "weeks")
    ctDNA.plot$START_DATE_weeks <- as.double(ctDNA.plot$START_DATE)/7
    ctDNA.plot$START_DATE_weeks <- ctDNA.plot$START_DATE_weeks - 
      as.double (plot.dat[as.character(ctDNA.plot$PATIENT_ID),
                          "WEEKS_SINCE_DIAGNOSIS"])
    ctDNA.plot$Y.pos <- mp[as.numeric(factor(as.character(ctDNA.plot$PATIENT_ID),
                                             levels = row.names(plot.dat)[order(plot.dat[,"Max.time"], decreasing = FALSE)]))]
    max.length = ddply (ctDNA.plot, ~PATIENT_ID, summarize, max.val = max(START_DATE_weeks, na.rm = TRUE))
    max.length$Y.pos= mp[as.numeric(factor(as.character(max.length$PATIENT_ID), 
                                           levels = row.names(plot.dat)[order(plot.dat[,"Max.time"], decreasing = FALSE)]))]
    
  }else{
    load ("/Volumes/Samwise/projects/INSPIRE/Review/ctDNA/NeteraMeeting_07052019/ESMO_siwmmer_v2_06082019.Rdata")
    timeline$colors[timeline$colors=="#FFFF99"] <- "#FDD901"
  }
  
  # subset patients if given, if not, show entire dataset
  if (!is.na(pt.select)){
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
  
  # plotting start here
  # Plot axis labels
  x.lab <- "Time since initiation of therapy, wks"
  y.lab <- "Subjects recieved study drug, pts"
  x.lim <- c(0,round_any(as.numeric( max(plot.dat$Max.time, na.rm = TRUE)+10), 20, f = ceiling))
  print(x.lim)
  opar <- par(lwd = 0.3)
  
  # calculate variables that are re-used by multiple functions
  n.bar <- nrow(plot.dat)
  bar.gap <- 0.2
  bar.width <- 1
  y.max <- n.bar*bar.width + (1+n.bar)*bar.gap
  
  #plot each bar
  par(mar = c(6.15 ,5, 5.2, 2.1))

  plot(NA, 
       ylim = c(0,y.max), 
       xlim = x.lim, 
       xlab = "", 
       ylab = "",
       yaxs="i",
       xaxt = "n", 
       yaxt = "n")
  
  segments(x0 = timeline$interval.start, 
           x1 = timeline$interval.end, 
           y0 = mp[timeline$Y.pos],
           col = timeline$colors,
           lwd = 5)
  # plot arrow for those continuing treatment
  points (
    x = plot.dat[plot.dat[,"Still on treatment (y or n)"]=="Y","Max.time"]+1.5,
    y = mp[plot.dat[,"Still on treatment (y or n)"]=="Y"],
    pch = -c(8594),
    cex = 1,
    col = "black")
  
  # plot symbol for those off trial not due to progression
  if(decision){
    #reasons.keep = c("Death", "Intercurrent Illness", "Investigator's Decision",
    #                 "Toxicity", "Withdrew consent", "Withdrew Consent")
    reasons.keep <- "Death"
    points (
      x = plot.dat[plot.dat[,"REASON_OFF_TRIAL"] %in% reasons.keep,"Weeks.on.treatment"],
      y = mp[plot.dat[,"REASON_OFF_TRIAL"] %in% reasons.keep],
      pch = 4,
      lwd = 2,
      cex = 0.70,
      col = "black")
  }
  
  # plot deaths
  points (
    x = plot.dat[as.character(plot.dat[,"OS_STATUS"]) %in% "DECEASED","Max.time"],
    y = mp[as.character(plot.dat[,"OS_STATUS"]) %in% "DECEASED"],
    pch = 4,
    lwd = 2,
    cex = 0.70,
    col = "black")                                                    
  
  points(
    x = ifelse(ctDNA.plot$START_DATE_weeks<0, 0, ctDNA.plot$START_DATE_weeks),
    y = ctDNA.plot$Y.pos,
    col = "black",
    bg = as.character(ctDNA.plot$DetectedStatus),
    pch = 21,
    cex = 0.70
  )
  
  title(xlab = x.lab, main = "RECIST 1.1 Clincal Response", cex.lab = 0.7, line = 2.50, cex.main = 0.7)
  axis (side = 1, lwd = 1, labels = T, lwd.tick=1, xaxs="i", cex.axis = 0.7)
  axis (side = 2, lwd = 1, labels = row.names(plot.dat), at = mp, lwd.tick=1, xaxs="i", cex.axis = 0.5, las = 2)
  box(lwd = 0.75)
  
  swimmer.legend = list(
    labels = c("Stable disease",
               "Partial response",
               "Complete response",
               "Disease progression",
               "No response data",
               "Continued response",
               "Death"),
    ltys = c(1,1,1,1,1, NA, NA),
    lwds = c(2.5,2.5,2.5,2.5,2.5,NA,1.5),
    pchs = c(NA,NA,NA,NA,NA,-8594,4),
    pt.cexs = c(NA,NA,NA,NA,NA,1.5,1.3),
    pt.cols = c(barcol,"black","black")
  );
  
  #plot legend
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
  
  return (list(plot.dat = plot.dat, tm.all.tps.df = tm.all.tps.df, mp=mp,covariates=covariates,
               ctDNAsummary = ctDNAsummary,
               n.bar =n.bar, bar.gap = bar.gap, bar.width = bar.width, y.max = y.max, y.lab = y.lab,
               annotcol = as.character(factor(plot.dat[,"INSPIRE COHORT"], levels = names(cancer.cols), labels = cancer.cols))))
}