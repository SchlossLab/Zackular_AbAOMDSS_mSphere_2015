source("code/rf_baseline_analysis.R")

shared_file <- "data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared"
rel_abund <- get_rel_abund(paste0(PROJHOME, "/", shared_file), "DSS")
rel_abund <- rel_abund[!grepl("dneg", rownames(rel_abund)),]

del_met <- "VancStrep"
del_met_rabund <- rel_abund[grep(del_met, rownames(rel_abund)),]

del_vanc <- "StrepMetro"
del_vanc_rabund <- rel_abund[grep(del_vanc, rownames(rel_abund)),]

control <- "NoAbs"
control_rabund <- rel_abund[grep(control, rownames(rel_abund)),]

del_rabund <- rbind(del_met_rabund, del_vanc_rabund, control_rabund)

del_days <- sapply(rownames(del_rabund), get_day)
del_days[del_days==70] <- 73
del_days[del_days==68] <- 69
del_days[del_days==45] <- 46
del_days[del_days==39] <- 41
del_treat <- sapply(rownames(del_rabund), get_treatment)

samples <- c(0, 21, 25, 41, 46, 69, 73)

good <- del_days %in% samples
del_days <- del_days[good]
del_treat <- del_treat[good]
del_rabund <- del_rabund[good,]

otu03 <- aggregate(del_rabund$Otu00003, by=list(del_days, del_treat), summary)
del_met_otu03 <- otu03[otu03$Group.2 == del_met,]
del_vanc_otu03 <- otu03[otu03$Group.2 == del_vanc,]
control_otu03 <- otu03[otu03$Group.2 == control,]

otu06 <- aggregate(del_rabund$Otu00006, by=list(del_days, del_treat), summary)
del_met_otu06 <- otu06[otu06$Group.2 == del_met,]
del_vanc_otu06 <- otu06[otu06$Group.2 == del_vanc,]
control_otu06 <- otu06[otu06$Group.2 == control,]

otu19 <- aggregate(del_rabund$Otu00019, by=list(del_days, del_treat), summary)
del_met_otu19 <- otu19[otu19$Group.2 == del_met,]
del_vanc_otu19 <- otu19[otu19$Group.2 == del_vanc,]
control_otu19 <- otu19[otu19$Group.2 == control,]

jitter <- 0.5

#par(mar=c(0.5,3,0.5,0.5))
plot(NA, ylim=c(0, 0.75), xlim=c(0,75), xlab="", ylab="", axes=F)
polygon(c(5,5,10,10),c(80,-10,-10,80), col='light blue', border=NA)   # DSS 1 box
polygon(c(26,26,31,31),c(80,-10,-10,80), col='light blue', border=NA) # DSS 2 box
polygon(c(47,47,52,52),c(80,-10,-10,80), col='light blue', border=NA) # DSS 3 box

segments(x0=del_met_otu03$Group.1+jitter, y0=del_met_otu03$x[,"Min."],
            y1=del_met_otu03$x[,"Max."], lwd=2, col=clrs[del_met_otu03$Group.2])
points(del_met_otu03$Group.1+jitter, del_met_otu03$x[,"Median"], col=clrs[del_met_otu03$Group.2], type="l")
points(del_met_otu03$Group.1+jitter, del_met_otu03$x[,"Median"], pch=21, col=clrs[del_met_otu03$Group.2], bg="white", cex=1.5)

segments(x0=del_vanc_otu03$Group.1-jitter, y0=del_vanc_otu03$x[,"Min."],
            y1=del_vanc_otu03$x[,"Max."], lwd=2, col=clrs[del_vanc_otu03$Group.2])
points(del_vanc_otu03$Group.1-jitter, del_vanc_otu03$x[,"Median"], col=clrs[del_vanc_otu03$Group.2], type="l")
points(del_vanc_otu03$Group.1-jitter, del_vanc_otu03$x[,"Median"], pch=21, col=clrs[del_vanc_otu03$Group.2], bg="white", cex=1.5)

segments(x0=control_otu03$Group.1, y0=control_otu03$x[,"Min."],
            y1=control_otu03$x[,"Max."], lwd=2, col=clrs[control_otu03$Group.2])
points(control_otu03$Group.1, control_otu03$x[,"Median"], col=clrs[control_otu03$Group.2], type="l")
points(control_otu03$Group.1, control_otu03$x[,"Median"], pch=21, col=clrs[control_otu03$Group.2], bg="white", cex=1.5)

axis(2, at=seq(0,0.75,0.25), label=seq(0,75,25), las=2)
text(x=78, y=0.73, label="Enterobacteriaceae (OTU 3)", pos=2, font=2)
box()


plot(NA, ylim=c(0, 0.75), xlim=c(0,75), xlab="", ylab="", axes=F)
polygon(c(5,5,10,10),c(80,-10,-10,80), col='light blue', border=NA)   # DSS 1 box
polygon(c(26,26,31,31),c(80,-10,-10,80), col='light blue', border=NA) # DSS 2 box
polygon(c(47,47,52,52),c(80,-10,-10,80), col='light blue', border=NA) # DSS 3 box

segments(x0=del_met_otu06$Group.1+jitter, y0=del_met_otu06$x[,"Min."],
            y1=del_met_otu06$x[,"Max."], lwd=2, col=clrs[del_met_otu06$Group.2])
points(del_met_otu06$Group.1+jitter, del_met_otu06$x[,"Median"], col=clrs[del_met_otu06$Group.2], type="l")
points(del_met_otu06$Group.1+jitter, del_met_otu06$x[,"Median"], pch=21, col=clrs[del_met_otu06$Group.2], bg="white", cex=1.5)

segments(x0=del_vanc_otu06$Group.1-jitter, y0=del_vanc_otu06$x[,"Min."],
            y1=del_vanc_otu06$x[,"Max."], lwd=2, col=clrs[del_vanc_otu06$Group.2])
points(del_vanc_otu06$Group.1-jitter, del_vanc_otu06$x[,"Median"], col=clrs[del_vanc_otu06$Group.2], type="l")
points(del_vanc_otu06$Group.1-jitter, del_vanc_otu06$x[,"Median"], pch=21, col=clrs[del_vanc_otu06$Group.2], bg="white", cex=1.5)

segments(x0=control_otu06$Group.1, y0=control_otu06$x[,"Min."],
            y1=control_otu06$x[,"Max."], lwd=2, col=clrs[control_otu06$Group.2])
points(control_otu06$Group.1, control_otu06$x[,"Median"], col=clrs[control_otu06$Group.2], type="l")
points(control_otu06$Group.1, control_otu06$x[,"Median"], pch=21, col=clrs[control_otu06$Group.2], bg="white", cex=1.5)

axis(2, at=seq(0,0.75,0.25), label=seq(0,75,25), las=2)
text(x=78, y=0.73, label="Lactobacillus (OTU 6)", pos=2, font=2)
box()


plot(NA, ylim=c(0, 0.25), xlim=c(0,75), xlab="", ylab="", axes=F)
polygon(c(5,5,10,10),c(80,-10,-10,80), col='light blue', border=NA)   # DSS 1 box
polygon(c(26,26,31,31),c(80,-10,-10,80), col='light blue', border=NA) # DSS 2 box
polygon(c(47,47,52,52),c(80,-10,-10,80), col='light blue', border=NA) # DSS 3 box

segments(x0=del_met_otu19$Group.1+jitter, y0=del_met_otu19$x[,"Min."],
            y1=del_met_otu19$x[,"Max."], lwd=2, col=clrs[del_met_otu19$Group.2])
points(del_met_otu19$Group.1+jitter, del_met_otu19$x[,"Median"], col=clrs[del_met_otu19$Group.2], type="l")
points(del_met_otu19$Group.1+jitter, del_met_otu19$x[,"Median"], pch=21, col=clrs[del_met_otu19$Group.2], bg="white", cex=1.5)

segments(x0=del_vanc_otu19$Group.1-jitter, y0=del_vanc_otu19$x[,"Min."],
            y1=del_vanc_otu19$x[,"Max."], lwd=2, col=clrs[del_vanc_otu19$Group.2])
points(del_vanc_otu19$Group.1-jitter, del_vanc_otu19$x[,"Median"], col=clrs[del_vanc_otu19$Group.2], type="l")
points(del_vanc_otu19$Group.1-jitter, del_vanc_otu19$x[,"Median"], pch=21, col=clrs[del_vanc_otu19$Group.2], bg="white", cex=1.5)

segments(x0=control_otu19$Group.1, y0=control_otu19$x[,"Min."],
            y1=control_otu19$x[,"Max."], lwd=2, col=clrs[control_otu19$Group.2])
points(control_otu19$Group.1, control_otu19$x[,"Median"], col=clrs[control_otu19$Group.2], type="l")
points(control_otu19$Group.1, control_otu19$x[,"Median"], pch=21, col=clrs[control_otu19$Group.2], bg="white", cex=1.5)

axis(2, at=seq(0,0.25,0.05), label=seq(0,25,5), las=2)
text(x=78, y=0.24, label="Prevotella (OTU 19)", pos=2, font=2)
box()

axis(1, at=seq(0,75,25), label=seq(0,75,25))
