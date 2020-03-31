# R code for Figures 1, 3, 5A, supplementary tables 2, 3, 4, 5, and
# supplementary figures 5, 7, 10

rm(list = ls())
library(Hmisc)
library(Epi)
library(reportRx)
library(forestplot)
library(survminer)
library(aod)
library(xtable)
library(pROC)

dat <- read.xlsx("Z:/Projects/INSPIRE/share/dat.xlsx", 1)
dim(dat)

# DNA1.cat: baseline ctDNA (dl), dichotomized (categorical variable)
# DNA1: baseline ctDNA (dl), as a continuous variable Change.cat:
# relative change in ctDNA between baseline and cycle 3, dichotomized
# Change: relative change in ctDNA, as a continuous variable DNA2:
# ctDNA measured at cycle 3 DNA2.cat: ctDNA measured at cycle 3,
# dichotomized PD_L1_percent_cat: PD-L1 percentage (categorical
# variable) log.TMB: logarithm of TMB

dat$change.cat <- relevel(dat$change.cat, ref = c("Increase from baseline"))

# pts with baseline ctDNA
allpts_withbaseline <- subset(dat, !is.na(dat$DNA1))
dim(allpts_withbaseline)

# pts with both baseline and C3 ctDNA
pts_withbaseline_C3_74 <- subset(dat, !is.na(dat$DNA1) & !is.na(dat$DNA2))
dim(pts_withbaseline_C3_74)

# 2 pts with baseline ctDNA=0; their delta ctDNA cannot be calculated
pts_withbaseline_C3 <- subset(dat, !is.na(dat$DNA1) & !is.na(dat$change))
dim(pts_withbaseline_C3)

getnumatrisk <- function(times, timevar) {
  atrisk <- NULL
  for (t in times) atrisk <- c(atrisk, sum(timevar >= t))
  return(atrisk)
}

# Figure 1 KM and Forest plot

mycol <- c("black", "red", "blue")
mylty <- c(1, 2, 4)

plot.km <- function(years = 36, by = 6, data = allpts_withbaseline, time = "OSTIME_Months", 
                    response = "OSevent", group = "DNA1.cat", lab = "Overall Survival (%)", 
                    cohort = "cohort", legend.position = "topright", xlab = "Months") {
  fit <- survfit(as.formula(paste("Surv(", time, ",", response, ")~", 
                                  paste(group, collapse = "+"), sep = "")), data = data)
  m <- coxph(as.formula(paste("Surv(", time, ",", response, ")~", 
                              paste(group, collapse = "+"), " + ", paste(cohort, collapse = "+"), 
                              sep = "")), data = data)
  HR <- round(coefficients(summary(m))[1, 2], digits = 2)
  pval <- Hmisc::format.pval(coef(summary(m))[1, 5], digits = 3, eps = 0.001)
  
  cats <- sort(unique(data[, group]))
  
  par(mfrow = c(1, 1), mar = c(7, 5, 1, 1))
  grid <- seq(0, years, by = by)
  
  plot(fit, xlab = "", ylab = "", lwd = 2, xlim = c(0, years), yaxt = "n", 
       xaxt = "n", mark.time = T, col = mycol[1:length(cats)], lty = 1:length(cats))
  axis(2, at = seq(0, 1, 0.2), las = TRUE, lab = paste(seq(0, 1, 0.2) * 
                                                         100, sep = ""))
  axis(1, at = grid, las = TRUE)
  mtext(xlab, side = 1, line = 2.5)
  mtext(lab, side = 2, line = 3)
  mtext("N at risk:", side = 1, line = 4, at = -4)
  for (i in 1:length(cats)) mtext(getnumatrisk(grid, data[, time][data[, group] == cats[i]]), 
                                  side = 1, line = 3 + i, at = grid, col = mycol[i])
  
  legend(legend.position, col = mycol[1:length(cats)], lwd = 2, legend = cats, 
         lty = 1:length(cats), bty = "n")
}

# 1A OS - baseline ctDNA
temp <- allpts_withbaseline

plot.km(years = 36, by = 6, data = temp, time = "OSTIME_Months", response = "OSevent", 
        group = "DNA1.cat", lab = "Overall Survival (%)")

# 1B PFS - baseline ctDNA
plot.km(years = 36, by = 6, data = temp, time = "PFSTIME_Months", response = "PFSevent", 
        group = "DNA1.cat", lab = "Progression-Free Survival (%)")

# 1C OS - baseline ctDNA by cohort
cohortlist <- c("A", "B", "C", "D", "E")
mat <- matrix(NA, ncol = 1, nrow = 5)
for (i in 1:5) {
  mat[i, ] <- mvsum(coxph(Surv(OSTIME_Months, OSevent) ~ DNA1.cat, 
                          data = subset(temp, temp$cohort == cohortlist[i])), 
                    data = subset(temp, temp$cohort == cohortlist[i]))[3, 2]
}

x6 <- mvsum(coxph(Surv(OSTIME_Months, OSevent) ~ DNA1.cat + cohort, data = temp), 
            data = temp)

tt <- rbind(mat, x6[3, 2])
tt1 <- gsub("\\(", ",", tt)
tt2 <- gsub("\\)", "", tt1)
tt3 <- strsplit(tt2, ",")
tt4 <- unlist(tt3)
tt5 <- as.numeric(tt4)
tt6 <- matrix(tt5, nrow = 3, ncol = 6, byrow = F)
tt6[, 4] <- NA
mean <- c(NA, tt6[1, ])
lower <- c(NA, tt6[2, ])
upper <- c(NA, tt6[3, ])
tt[4, ] <- "NA"

tabletext <- cbind(c("Cohort", "A", "B", "C", "D", "E", "Overall"), 
                   c("N", paste0(table(temp$cohort)), nrow(temp)), c("HR (95% CI)", tt))


#In cohort D, all the event times in the above-median group were earlier than those in the below-median group, and thus HR cannot be estimated
forestplot(labeltext = tabletext, mean = mean, lower = lower, upper = upper, 
           xlog = T, align = "c", is.summary = c(TRUE, rep(FALSE, 5), TRUE), 
           col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"), 
           vertices = TRUE, xlab = "HR (95% CI)", 
           grid = structure(c(1), gp = gpar(lty = 2, col = "red")), hrzl_lines = gpar(col = "#444444"), 
           xticks.digits = 1, txt_gp = fpTxtGp(label = gpar(fontfamily = ""), 
           xlab = gpar(cex = 1), ticks = gpar(cex = 1)))


# 1D PFS - baseline ctDNA by cohort
mat <- matrix(NA, ncol = 1, nrow = 5)
for (i in 1:5) {
  mat[i, ] <- mvsum(coxph(Surv(PFSTIME_Months, PFSevent) ~ DNA1.cat, 
                          data = subset(temp, temp$cohort == cohortlist[i])), 
                    data = subset(temp, temp$cohort == cohortlist[i]))[3, 2]
}

x6 <- mvsum(coxph(Surv(PFSTIME_Months, PFSevent) ~ DNA1.cat + cohort, data = temp), 
            data = temp)

tt <- rbind(mat, x6[3, 2])
tt1 <- gsub("\\(", ",", tt)
tt2 <- gsub("\\)", "", tt1)
tt3 <- strsplit(tt2, ",")
tt4 <- unlist(tt3)
tt5 <- as.numeric(tt4)
tt6 <- matrix(tt5, nrow = 3, ncol = 6, byrow = F)
tt6[, 4] <- NA
mean <- c(NA, tt6[1, ])
lower <- c(NA, tt6[2, ])
upper <- c(NA, tt6[3, ])
tt[4, ] <- "NA"

tabletext <- cbind(c("Cohort", "A", "B", "C", "D", "E", "Overall"), 
                   c("N", paste0(table(temp$cohort)), nrow(temp)), c("HR", tt))

#In cohort D, all the event times in the above-median group were earlier than those in the below-median group, and thus HR cannot be estimated

forestplot(labeltext = tabletext, mean = mean, lower = lower, upper = upper, 
           xlog = T, align = "c", is.summary = c(TRUE, rep(FALSE, 5), TRUE), 
           col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"), 
           vertices = TRUE, grid = structure(c(1),  gp = gpar(lty = 2, col = "red")), 
           xlab = "HR(95% CI)", hrzl_lines = gpar(col = "#444444"), 
           xticks.digits = 1, txt_gp = fpTxtGp(label = gpar(fontfamily = ""), 
           xlab = gpar(cex = 1), ticks = gpar(cex = 1)))



# Fig 3

# 3A OS - change in ctDNA
temp <- pts_withbaseline_C3

plot.km(years = 36, by = 6, data = temp, time = "OSTIME_Months_C3", response = "OSevent", 
        group = "change.cat", lab = "Overall Survival (%)", xlab = "Months from C3")

# 3B PFS - change in ctDNA Patient INS-B-002 and INS-D-008 had progressed before C3
plot.km(years = 36, by = 6, data = temp, time = "PFSTIME_Months_C3", response = "PFSevent", 
        group = "change.cat", lab = "Progression-Free Survival (%)", xlab = "Months from C3")

# 3C OS - change in ctDNA by cohort
mat <- matrix(NA, ncol = 1, nrow = 5)
for (i in 1:5) {
  mat[i, ] <- mvsum(coxph(Surv(OSTIME_Months_C3, OSevent) ~ change.cat, 
                          data = subset(temp, temp$cohort == cohortlist[i])), 
                    data = subset(temp, temp$cohort == cohortlist[i]))[3, 2]
}

x6 <- mvsum(coxph(Surv(OSTIME_Months_C3, OSevent) ~ change.cat + cohort, 
                  data = temp), data = temp)

tt <- rbind(mat, x6[3, 2])
tt1 <- gsub("\\(", ",", tt)
tt2 <- gsub("\\)", "", tt1)
tt3 <- strsplit(tt2, ",")
tt4 <- unlist(tt3)
tt5 <- as.numeric(tt4)
tt6 <- matrix(tt5, nrow = 3, ncol = 6, byrow = F)
mean <- c(NA, tt6[1, ])
lower <- c(NA, tt6[2, ])
upper <- c(NA, tt6[3, ])


tabletext <- cbind(c("Cohort", "A", "B", "C", "D", "E", "Overall"), 
                   c("N", paste0(table(temp$cohort)), nrow(temp)), c("HR (95% CI)", tt))

forestplot(labeltext = tabletext, mean = mean, lower = lower, upper = upper, 
           xlog = T, align = "c", is.summary = c(TRUE, rep(FALSE, 5), TRUE), 
           col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"), 
           vertices = TRUE, grid = structure(c(1), gp = gpar(lty = 2, col = "red")), 
           xlab = "HR(95% CI)", hrzl_lines = gpar(col = "#444444"), 
           xticks.digits = 1, txt_gp = fpTxtGp(label = gpar(fontfamily = ""), 
           xlab = gpar(cex = 1), ticks = gpar(cex = 1)))


# 3D PFS - change in ctDNA by cohort
mat <- matrix(NA, ncol = 1, nrow = 5)
for (i in 1:5) {
  mat[i, ] <- mvsum(coxph(Surv(PFSTIME_Months_C3, PFSevent) ~ change.cat, 
                          data = subset(temp, temp$cohort == cohortlist[i])), 
                    data = subset(temp, temp$cohort == cohortlist[i]))[3, 2]
}

x6 <- mvsum(coxph(Surv(PFSTIME_Months_C3, PFSevent) ~ change.cat + cohort, 
                  data = temp), data = temp)

tt <- rbind(mat, x6[3, 2])
tt1 <- gsub("\\(", ",", tt)
tt2 <- gsub("\\)", "", tt1)
tt3 <- strsplit(tt2, ",")
tt4 <- unlist(tt3)
tt5 <- as.numeric(tt4)
tt6 <- matrix(tt5, nrow = 3, ncol = 6, byrow = F)
mean <- c(NA, tt6[1, ])
lower <- c(NA, tt6[2, ])
upper <- c(NA, tt6[3, ])

temp1 <- subset(temp, temp$PFSTIME_Months_C3 >= 0)

tabletext <- cbind(c("Cohort", "A", "B", "C", "D", "E", "Overall"), 
                   c("N", paste0(table(temp1$cohort)), nrow(temp1)), c("HR", tt))


forestplot(labeltext = tabletext, mean = mean, lower = lower, upper = upper, 
           xlog = T, align = "c", is.summary = c(TRUE, rep(FALSE, 5), TRUE), 
           col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"), 
           vertices = TRUE, grid = structure(c(1), gp = gpar(lty = 2, col = "red")), 
           xlab = "HR(95% CI)", hrzl_lines = gpar(col = "#444444"), 
           xticks.digits = 1, txt_gp = fpTxtGp(label = gpar(fontfamily = ""), 
                                               xlab = gpar(cex = 1), ticks = gpar(cex = 1)))

# Fig 5A
pts_withbaseline_C3_73 <- subset(pts_withbaseline_C3, pts_withbaseline_C3$id != 
                                   "INS-D-012")

plot.km.noHR2 <- function(years = 36, by = 6, data = pts_withbaseline_C3_73, 
                          time = "OSTIME_Months_C3", response = "OSevent", group = "clearance.group", 
                          legend.x = 0, legend.y = 0.2) {
  fit <- survfit(as.formula(paste("Surv(", time, ",", response, ")~", 
                                  paste(group, collapse = "+"), sep = "")), data = data)
  cats <- c("With clearance", "No clearance, decrease from baseline", 
            "No clearance, increase from baseline")
  par(mfrow = c(1, 1), mar = c(7, 5, 1, 1))
  
  grid <- seq(0, years, by = by)
  
  plot(fit, xlab = "", ylab = "Overall Survival (%)", lwd = 2, xlim = c(0, years), 
       yaxt = "n", xaxt = "n", mark.time = T, col = c(1, 2, 4), lty = c(1, 2, 4), cex.lab = 1)
  axis(2, at = seq(0, 1, 0.2), las = TRUE, lab = paste(seq(0, 1, 0.2) * 
                                                         100, sep = ""))
  axis(1, at = grid, las = TRUE)
  mtext("Months on treatment", side = 1, line = 2.5, cex = 1)
  mtext("N at risk:", side = 1, line = 4, at = -4)
  for (i in 1:length(cats)) mtext(getnumatrisk(grid, data[, time][data[, group] == cats[i]]), 
                                  side = 1, line = 3 + i, at = grid)
  
  legend(legend.x, legend.y, col = c(4, 1, 2), lty = c(4, 1, 2), lwd = 2, 
         legend = c("With clearance", "No clearance, decrease from baseline", 
                    "No clearance, increase from baseline"), bty = "n")
}

plot.km.noHR2(years = 36, by = 6, data = pts_withbaseline_C3_73, time = "OSTIME_Months_C3", 
              response = "OSevent", group = "clearance.group")

# Supp Table 2. Patient demographics
pcovsum(allpts_withbaseline, covs = c("Age_at_diagnosis", "Sex", "Ethnicity1", 
                                      "cohort"), maincov = NULL, testcat = "Fisher", testcont = c("ANOVA"))

# Supp Table 3. UVA with baseline ctDNA

# OS - UVA
temp <- allpts_withbaseline
covs <- "DNA1.cat"
puvsum(c("OSTIME_Months", "OSevent"), covs = covs, data = temp, type = "coxph")

# PFS - UVA
puvsum(c("PFSTIME_Months", "PFSevent"), covs = covs, data = temp, type = "coxph")

# ORR - UVA
puvsum(c("Best_response1"), covs = covs, data = temp, type = "logistic")

# CBR - UVA
puvsum(c("CR_PR_SD_6_cycles"), covs = covs, data = temp, type = "logistic")


# Supp Table 4. UVA for change in ctDNA

# OS - UVA
temp <- pts_withbaseline_C3_74
covs <- c("change.cat")
puvsum(c("OSTIME_Months_C3", "OSevent"), covs = covs, data = temp, type = "coxph")

# PFS - UVA
puvsum(c("PFSTIME_Months_C3", "PFSevent"), covs = covs, data = temp, type = "coxph")

# ORR - UVA
puvsum(c("Best_response1"), covs = covs, data = temp, type = "logistic")

# CBR - UVA
puvsum(c("CR_PR_SD_6_cycles"), covs = covs, data = temp, type = "logistic")


# Supp Table 5. MVA with delta ctDNA

temp <- pts_withbaseline_C3_74

# OS - MVA
f2 <- coxph(Surv(OSTIME_Months_C3, OSevent) ~ change.cat + cohort + PD_L1_percent + 
              log.TMB, temp)
pmvsum(f2, temp)

# PFS - MVA
f2 <- coxph(Surv(PFSTIME_Months_C3, PFSevent) ~ change.cat + cohort + PD_L1_percent + 
              log.TMB, temp)
pmvsum(f2, temp)

# ORR - MVA
cat("adjusting for PD-L1")
f <- glm(Best_response1 ~ change.cat + PD_L1_percent, data = temp, family = "binomial")
pmvsum(f, temp)

cat("adjusting for TMB")
f <- glm(Best_response1 ~ change.cat + log.TMB, data = temp, family = "binomial")
pmvsum(f, temp)

# CBR - MVA
cat("adjusting for PD-L1")
f <- glm(CR_PR_SD_6_cycles ~ change.cat + PD_L1_percent, data = temp, family = "binomial")
pmvsum(f, temp)

cat("adjusting for TMB")
f <- glm(CR_PR_SD_6_cycles ~ change.cat + log.TMB, data = temp, family = "binomial")
pmvsum(f, temp)


# Supp Fig 5 baseline ctDNA and target lesion size at baseline

allpts_withbaseline$DNA1.original <- allpts_withbaseline$DNA1 * 100
x <- allpts_withbaseline$DNA1.original
y <- allpts_withbaseline$MALIGN_SUM_DIAM_VOL
plot(x, y, xlab = "Mean Tumour Molecules per mL Plasma", ylab = "Target lesion size at baseline (mm)")
p <- cor.test(x, y, method = c("spearman"), exact = F)$p.value
c <- cor.test(x, y, method = c("spearman"), exact = F)$estimate
legend("topright", legend = paste0("Spearman correlation = ", round(c, dig = 2), "; p = ", 
                                   round(p, digits = 2)))


# Supp Fig 7 scatter plots
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))

x <- allpts_withbaseline$DNA1.original
y <- allpts_withbaseline$log.TMB
plot(x, y, xlab = "Mean Tumour Molecules per mL Plasma", ylab = "log(TMB)", 
     yaxt = "n", ylim = c(-1, 5.5))
axis(side = 2, at = c(-1, 0:5), labels = c(-1, 0:5))
p <- cor.test(x, y, method = c("spearman"), exact = F)$p.value
c <- cor.test(x, y, method = c("spearman"), exact = F)$estimate
legend("topright", legend = paste0("Spearman correlation = ", round(c, dig = 2), "; p = ", 
                                   round(p, digits = 2)), bg = "transparent")

x <- allpts_withbaseline$DNA1.original
y <- allpts_withbaseline$PD_L1_percent
plot(x, y, xlab = "Mean Tumour Molecules per mL Plasma", ylab = "PD-L1 percentage (%)", 
     ylim = c(0, 105))
p <- cor.test(x, y, method = c("spearman"), exact = F)$p.value
c <- cor.test(x, y, method = c("spearman"), exact = F)$estimate
legend("topright", legend = paste0("Spearman correlation = ", round(c, dig = 2), "; p = ", 
                                   round(p, digits = 2)), bg = "transparent")

sub11 <- pts_withbaseline_C3
xnew <- ifelse(sub11$change * 100 < 0, sub11$change * 100, sub11$change * 
                 100/10)
plot(xnew, sub11$log.TMB, ylab = "log(TMB)", xlab = "% Change in ctDNA", 
     yaxt = "n", xaxt = "n", ylim = c(-1, 6), xlim = c(-150, 150))
axis(side = 2, at = c(-1, 0:7), labels = c(-1, 0:7))
axis(side = 1, at = seq(-150, 150, by = 50), labels = c(-150, -100, -50, 
                                                        0, 500, 1000, 1500))
p <- cor.test(sub11$log.TMB, sub11$change * 100, method = c("spearman"), 
              exact = F)$p.value
c <- cor.test(sub11$log.TMB, sub11$change * 100, method = c("spearman"), 
              exact = F)$estimate
legend("topright", legend = paste0("Spearman correlation = ", round(c, dig = 2), "; p = ", 
                                   round(p, digits = 2)), bg = "transparent")


plot(xnew, sub11$PD_L1_percent, ylab = "PD-L1 percentage (%)", xlab = "% Change in ctDNA", 
     ylim = c(0, 115), xlim = c(-150, 150), xaxt = "n")
axis(side = 1, at = seq(-150, 150, by = 50), labels = c(-150, -100, -50, 
                                                        0, 500, 1000, 1500))
p <- cor.test(sub11$change * 100, sub11$PD_L1_percent, method = c("spearman"), 
              exact = F)$p.value
c <- cor.test(sub11$change * 100, sub11$PD_L1_percent, method = c("spearman"), 
              exact = F)$estimate
legend("topright", legend = paste0("Spearman correlation = ", round(c, dig = 2), "; p = ", 
                                   round(p, digits = 2)), bg = "transparent")

dev.off()

# Supp Fig 10 Landmark analysis

# Clearance is assessed at Cycles 3, 6, 9 and 12.
mycol <- c("blue", "black", "red")
mylty <- c(4, 1, 2)

plot.km.noHR <- function(years = 36, by = 6, data = allpts_withbaseline, 
                         time = "OSTIME_Months", response = "OSevent", group = "clearance.group", 
                         lab = "Overall Survival (%)", xlab = "Months from C3 landmark") {
  fit <- survfit(as.formula(paste("Surv(", time, ",", response, ")~", 
                                  paste(group, collapse = "+"), sep = "")), data = data)
  cats <- sort(unique(data[, group]))
  
  grid <- seq(0, years, by = by)
  
  plot(fit, xlab = "", ylab = "", lwd = 2, xlim = c(0, years), yaxt = "n", 
       xaxt = "n", mark.time = T, col = mycol[1:length(cats)], lty = mylty[1:length(cats)])
  axis(2, at = seq(0, 1, 0.2), las = TRUE, lab = paste(seq(0, 1, 0.2) * 
                                                         100, sep = ""))
  axis(1, at = grid, las = TRUE)
  mtext(xlab, side = 1, line = 2.5)
  mtext(lab, side = 2, line = 3)
  mtext("N at risk:", side = 1, line = 4, at = -8)
  for (i in 1:length(cats)) mtext(getnumatrisk(grid, data[, time][data[, group] == cats[i]]), 
                                  side = 1, line = 3 + i, at = grid, col = mycol[i])
  
}

par(oma = c(4, 1, 1, 1), mfrow = c(2, 2), mar = c(7, 5, 0.5, 0.5))

temp <- subset(allpts_withbaseline, !is.na(allpts_withbaseline$DNA2))
temp$group <- ifelse(temp$DNA2 == 0, "Clearance", "N/A")
temp$group <- ifelse(temp$group == "N/A" & temp$DNA2 > temp$DNA1, "Increase from baseline", 
                     as.character(temp$group))
temp$group <- ifelse(temp$group == "N/A" & temp$DNA2 < temp$DNA1, "Decrease from baseline", 
                     as.character(temp$group))
temp$group <- as.factor(temp$group)
plot.km.noHR(years = 36, by = 6, data = temp, time = "OSTIME_Months_C3", 
             response = "OSevent", group = "group")

temp <- subset(allpts_withbaseline, !is.na(allpts_withbaseline$C6_ctDNA))
temp$group <- ifelse(temp$C6_ctDNA == 0, "Clearance", "N/A")
temp$group <- ifelse(temp$group == "N/A" & temp$C6_ctDNA > temp$DNA1, "Increase from baseline", 
                     as.character(temp$group))
temp$group <- ifelse(temp$group == "N/A" & temp$C6_ctDNA < temp$DNA1, "Decrease from baseline", 
                     as.character(temp$group))
temp$group <- as.factor(temp$group)
plot.km.noHR(years = 36, by = 6, data = temp, time = "OSTIME_Months_C6", 
             response = "OSevent", group = "group", xlab = "Months from C6 landmark")


temp <- subset(allpts_withbaseline, !is.na(allpts_withbaseline$C9_ctDNA))
temp$group <- ifelse(temp$C9_ctDNA == 0, "Clearance", "N/A")
temp$group <- ifelse(temp$group == "N/A" & temp$C9_ctDNA > temp$DNA1, "Increase from baseline", 
                     as.character(temp$group))
temp$group <- ifelse(temp$group == "N/A" & temp$C9_ctDNA < temp$DNA1, "Decrease from baseline", 
                     as.character(temp$group))
temp$group <- as.factor(temp$group)
plot.km.noHR(years = 36, by = 6, data = temp, time = "OSTIME_Months_C9", 
             response = "OSevent", group = "group", xlab = "Months from C9 landmark")

temp <- subset(allpts_withbaseline, !is.na(allpts_withbaseline$C12_ctDNA))
temp$group <- ifelse(temp$C12_ctDNA == 0, "Clearance", "N/A")
temp$group <- ifelse(temp$group == "N/A" & temp$C12_ctDNA > temp$DNA1, 
                     "Increase from baseline", as.character(temp$group))
temp$group <- ifelse(temp$group == "N/A" & temp$C12_ctDNA < temp$DNA1, 
                     "Decrease from baseline", as.character(temp$group))
temp$group <- as.factor(temp$group)
plot.km.noHR(years = 36, by = 6, data = temp, time = "OSTIME_Months_C12", 
             response = "OSevent", group = "group", xlab = "Months from C12 landmark")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", xpd = TRUE, col = mycol, lty = mylty, lwd = 2, 
       legend = c("Clearance", "Decrease from baseline", "Increase from baseline"), horiz = T, bty = "n")

dev.off()
