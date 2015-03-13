#cairo_pdf("test.pdf", width=4, height=7)

source("code/rf_baseline_analysis.R")

tax_file <- "data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.5.cons.taxonomy"
tax <- read.table(tax_file, header=T, row.names=1)
tax$Taxonomy <- gsub("^\\w*\\(\\d*\\);", "", tax$Taxonomy)
tax$Taxonomy <- gsub("^([^\\(]*)\\(.*", "\\1", tax$Taxonomy)

shared_file <- "data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.shared"
shared <- read_shared_file(shared_file, "DSS")

sample_names <- rownames(shared)
treatment <- sapply(sample_names, get_treatment)
day <- sapply(sample_names, get_day)

baseline_shared <- shared[day==0,]
baseline_treatment <- treatment[day==0]

baseline_rabund <- 100 * prop.table(as.matrix(baseline_shared), 1)
ordered_otus <- order(apply(baseline_rabund, 2, median), decreasing=T)

ordered_treatments <- c("NoAbs", "Vanc", "Strep", "VancStrep", "Metro",
                        "VancMetro", "AllAbs", "StrepMetro")


otu_plot <- function(rabunds, ordered_treatments, otu_label){

    means <- aggregate(rabunds, by=list(baseline_treatment), mean)$x
    o <- c(3, 6, 4, 8, 2, 7, 1, 5)

    par(mar=c(0.75, 2, 0.75, 0.5))

    plot(NA, xlim=c(0.75,8.25), xlab="", ylim=c(0,max(rabunds)), axes=F, cex.lab=1.4)

    segments(x0=seq(1:8)-0.4, x1=seq(1:8)+0.4, y0=means[o], y1=means[o], lwd=3)

    box()
    axis(2, las=2)

    for(i in 1:length(ordered_treatments)){
        stripchart(at=i, x=rabunds[baseline_treatment==ordered_treatments[i]],
        vertical=T, method="jitter", col=clrs[ordered_treatments[i]],
        pch=pch[ordered_treatments[i]], jitter=0.25, cex=1.5, add=T)
    }

    text(x=0.35, y=1.12*max(rabunds), label=otu_label, pos=4, font=2, cex=1, xpd=TRUE)

}

design <- matrix(c(1,2,3,4), ncol=1)
design <- cbind(c(5,5,5,5), design)
design <- rbind(design, c(0,6))

layout(design, width=c(0.13, 1), height=c(1,1,1,1,1.2))

otu_plot(baseline_rabund[,ordered_otus[1]], ordered_treatments, tax[ordered_otus[1],"Taxonomy"])
otu_plot(baseline_rabund[,ordered_otus[2]], ordered_treatments, tax[ordered_otus[2],"Taxonomy"])
otu_plot(baseline_rabund[,ordered_otus[3]], ordered_treatments, tax[ordered_otus[3],"Taxonomy"])
otu_plot(baseline_rabund[,ordered_otus[4]], ordered_treatments, tax[ordered_otus[4],"Taxonomy"])

plot.new()
text(x=0.7, y=0.5, label="Relative abundance (%)", cex=1.5, srt=90)
text(label="B", x=0.0, y=1.03, font=2.5, cex=2.5, xpd=TRUE)

plot.new()
text(x=seq(0.05,1,length=8), y=1, labels[ordered_treatments], srt=60, cex=1.0, pos=2, xpd=TRUE)

#dev.off()
