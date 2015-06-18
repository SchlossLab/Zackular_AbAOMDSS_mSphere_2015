source("code/rf_baseline_analysis.R")

shared_file <- "data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared"
rel_abund <- get_rel_abund(paste0(PROJHOME, "/", shared_file), "DSS")
rel_abund <- rel_abund[!grepl("dneg", rownames(rel_abund)),]

days <- sapply(rownames(rel_abund), get_day)
days[days==70] <- 73
days[days==68] <- 69
days[days==45] <- 46
days[days==39] <- 41

treat <- sapply(rownames(rel_abund), get_treatment)

samples <- c(0, 21, 25, 41, 46, 69, 73)

good <- days %in% samples
days <- days[good]
treat <- treat[good]
rel_abund <- rel_abund[good,]

get_top_delta <- function(treatment, n=5){

    treat_rel_abund <- rel_abund[treat==treatment,]
    treat_days <- days[treat==treatment]

    rabund_start_finish <- treat_rel_abund[treat_days==0 | treat_days==73,]
    days_start_finish <- treat_days[treat_days==0 | treat_days==73]
    et_start_finish <- gsub("_d\\d*", "", rownames(rabund_start_finish))

    intra_mouse_delta <- aggregate(rabund_start_finish,
                                    by=list(et_start_finish),
                                    function(x)max(x)-min(x))

    median_mouse_delta <- sort(apply(intra_mouse_delta[,-1], 2, median), decreasing=T)[1:n]

    names(median_mouse_delta)
}

all_otus <- get_top_delta("AllAbs")
metro_otus <- get_top_delta("Metro")
vanc_otus <- get_top_delta("VancStrep")
noabs_otus <- get_top_delta("NoAbs")
strepmetro_otus <- get_top_delta("StrepMetro")
vancmetro_otus <- get_top_delta("VancMetro")
vancstrep_otus <- get_top_delta("VancStrep")

composite_otus <- c(all_otus, metro_otus, vanc_otus, noabs_otus, strepmetro_otus,
                    vancmetro_otus, vancstrep_otus)
ranked_otus <- sort(table(composite_otus), decreasing=TRUE)




otu01 <- aggregate(rel_abund$Otu00001, by=list(days, treat), summary)
vancstrep_otu01 <- otu01[otu01$Group.2 == "VancStrep",]
strepmetro_otu01 <- otu01[otu01$Group.2 == "StrepMetro",]
noabx_otu01 <- otu01[otu01$Group.2 == "NoAbs",]

otu02 <- aggregate(rel_abund$Otu00002, by=list(days, treat), summary)
vancstrep_otu02 <- otu02[otu02$Group.2 == "VancStrep",]
strepmetro_otu02 <- otu02[otu02$Group.2 == "StrepMetro",]
noabx_otu02 <- otu02[otu02$Group.2 == "NoAbs",]

jitter <- 0.5

#par(mar=c(0.5,3,0.5,0.5))
plot(NA, ylim=c(0, 0.75), xlim=c(0,75), xlab="", ylab="", axes=F)
polygon(c(5,5,10,10),c(80,-10,-10,80), col='light blue', border=NA)   # DSS 1 box
polygon(c(26,26,31,31),c(80,-10,-10,80), col='light blue', border=NA) # DSS 2 box
polygon(c(47,47,52,52),c(80,-10,-10,80), col='light blue', border=NA) # DSS 3 box


segments(x0=vancstrep_otu01$Group.1+jitter, y0=vancstrep_otu01$x[,"Min."],
            y1=vancstrep_otu01$x[,"Max."], lwd=2, col=clrs[vancstrep_otu01$Group.2])
points(vancstrep_otu01$Group.1+jitter, vancstrep_otu01$x[,"Median"],
            col=clrs[vancstrep_otu01$Group.2], type="l")
points(vancstrep_otu01$Group.1+jitter, vancstrep_otu01$x[,"Median"],
            pch=21, col=clrs[vancstrep_otu01$Group.2], bg="white", cex=1.5)

segments(x0=strepmetro_otu01$Group.1-jitter, y0=strepmetro_otu01$x[,"Min."],
            y1=strepmetro_otu01$x[,"Max."], lwd=2, col=clrs[strepmetro_otu01$Group.2])
points(strepmetro_otu01$Group.1-jitter, strepmetro_otu01$x[,"Median"],
            col=clrs[strepmetro_otu01$Group.2], type="l")
points(strepmetro_otu01$Group.1-jitter, strepmetro_otu01$x[,"Median"], pch=21,
            col=clrs[strepmetro_otu01$Group.2], bg="white", cex=1.5)

segments(x0=noabx_otu01$Group.1-jitter, y0=noabx_otu01$x[,"Min."],
            y1=noabx_otu01$x[,"Max."], lwd=2, col=clrs[noabx_otu01$Group.2])
points(noabx_otu01$Group.1-jitter, noabx_otu01$x[,"Median"],
            col=clrs[noabx_otu01$Group.2], type="l")
points(noabx_otu01$Group.1-jitter, noabx_otu01$x[,"Median"], pch=21,
            col=clrs[noabx_otu01$Group.2], bg="white", cex=1.5)

axis(2, at=seq(0,0.75,0.25), label=seq(0,75,25), las=2)
text(x=78, y=0.73, label="Lactobacillus (OTU 1)", pos=2, font=2)
box()





plot(NA, ylim=c(0, 0.75), xlim=c(0,75), xlab="", ylab="", axes=F)
polygon(c(5,5,10,10),c(80,-10,-10,80), col='light blue', border=NA)   # DSS 1 box
polygon(c(26,26,31,31),c(80,-10,-10,80), col='light blue', border=NA) # DSS 2 box
polygon(c(47,47,52,52),c(80,-10,-10,80), col='light blue', border=NA) # DSS 3 box

segments(x0=vancstrep_otu02$Group.1+jitter, y0=vancstrep_otu02$x[,"Min."],
            y1=vancstrep_otu02$x[,"Max."], lwd=2, col=clrs[vancstrep_otu02$Group.2])
points(vancstrep_otu02$Group.1+jitter, vancstrep_otu02$x[,"Median"],
            col=clrs[vancstrep_otu02$Group.2], type="l")
points(vancstrep_otu02$Group.1+jitter, vancstrep_otu02$x[,"Median"], pch=21,
            col=clrs[vancstrep_otu02$Group.2], bg="white", cex=1.5)

segments(x0=strepmetro_otu02$Group.1-jitter, y0=strepmetro_otu02$x[,"Min."],
            y1=strepmetro_otu02$x[,"Max."], lwd=2, col=clrs[strepmetro_otu02$Group.2])
points(strepmetro_otu02$Group.1-jitter, strepmetro_otu02$x[,"Median"],
            col=clrs[strepmetro_otu02$Group.2], type="l")
points(strepmetro_otu02$Group.1-jitter, strepmetro_otu02$x[,"Median"], pch=21,
            col=clrs[strepmetro_otu02$Group.2], bg="white", cex=1.5)

segments(x0=noabx_otu02$Group.1-jitter, y0=noabx_otu02$x[,"Min."],
            y1=noabx_otu02$x[,"Max."], lwd=2, col=clrs[noabx_otu02$Group.2])
points(noabx_otu02$Group.1-jitter, noabx_otu02$x[,"Median"],
            col=clrs[noabx_otu02$Group.2], type="l")
points(noabx_otu02$Group.1-jitter, noabx_otu02$x[,"Median"], pch=21,
            col=clrs[noabx_otu02$Group.2], bg="white", cex=1.5)

axis(2, at=seq(0,0.75,0.25), label=seq(0,75,25), las=2)
text(x=78, y=0.73, label="Enterobacteriaceae (OTU 2)", pos=2, font=2)
box()


axis(1, at=seq(0,75,25), label=seq(0,75,25))
