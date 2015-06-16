library("png")
source("code/rf_baseline_analysis.R")

cairopdf("results/figures/figure1.pdf", width=7.5, height=10.0)

design <- c(1,1,10,
            3,2,10,
            8,4,10,
            8,5,10,
            8,6,11,
            8,7,11,
            0,9,11)

design <- matrix(design, ncol=3, byrow=TRUE)

z <- layout(design, width=c(0.2, 1, 1.2), height=c(2,2,1,1,1,1,1.2))


source("code/aomdss_timeline.R")
text(label="A", x=-14, y=4.5, font=2, cex=2)


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




tumors <- read.table(file="data/process/tumor_counts.tsv", header=T, row.names=1)

#get the mean number of tumors in each treatment
summarize_counts <- aggregate(tumors$n_tumors, by=list(baseline_treatment), mean)
means <- summarize_counts$x
names(means) <- summarize_counts$Group.1


#order the treatments by the number of tumors
decreasing_order <- order(means, decreasing=T)
ordered_treatments <- levels(as.factor(baseline_treatment))[decreasing_order]

counts <- table(baseline_treatment)

labels <- c(NoAbs = paste0("No antibiotics (N=",counts["NoAbs"],")"),
    AllAbs = paste0("All antibiotics (N=",counts["AllAbs"] , ")"),
    Metro = paste0("Metronidazole (N=",counts["Metro"] , ")"),
    VancStrep = paste0("\U0394Metronidazole (N=", counts["Metro"], ")"),
    Strep = paste0("Streptomycin (N=",counts["Strep"] , ")"),
    VancMetro = paste0("\U0394Streptomycin (N=",counts["Strep"] , ")"),
    Vanc = paste0("Vancomycin", " (N=",counts["Vanc"] , ")"),
    StrepMetro = paste0("\U0394Vancomycin (N=",counts["StrepMetro"], ")") )

par(mar=c(0.75, 2, 0.75, 0.5))
plot(NA, xlim=c(0.75,8.25), xlab="", ylim=c(0,25), axes=F, cex.lab=1.4)
segments(x0=seq(1:8)-0.4, x1=seq(1:8)+0.4, y0=means[decreasing_order], y1=means[decreasing_order], lwd=3)
box()
axis(2, las=2)
#axis(1, at=1:8, label=labels[ordered_treatments], las=2, tick=F)

#build a stripchart
for(i in 1:length(ordered_treatments)){
    stripchart(at=i, x=tumors[baseline_treatment==ordered_treatments[i], "n_tumors"],
    vertical=T, method="jitter", col=clrs[ordered_treatments[i]],
    pch=pch[ordered_treatments[i]], jitter=0.25, cex=1.0, add=T)
}
segments(x0=1,x1=8, y0=25, y1=25, lwd=2, col="gray")
segments(x0=1,x1=1, y0=25, y1=24, lwd=2)
text(x=c(3,4,5,6,7,8), y=rep(25, 6), label="*", cex=1.5)

segments(x0=2,x1=8, y0=22, y1=22, lwd=2, col="gray")
segments(x0=2,x1=2, y0=22, y1=21, lwd=2)
#

segments(x0=3,x1=8, y0=19, y1=19, lwd=2, col="gray")
segments(x0=3,x1=3, y0=19, y1=18, lwd=2)
text(x=c(7,8), y=rep(19, 2), label="*", cex=1.5)

segments(x0=4,x1=8, y0=16, y1=16, lwd=2, col="gray")
segments(x0=4,x1=4, y0=16, y1=15, lwd=2)
text(x=c(5,7,8), y=rep(16, 3), label="*", cex=1.5)

segments(x0=5,x1=8, y0=13, y1=13, lwd=2, col="gray")
segments(x0=5,x1=5, y0=13, y1=12, lwd=2)
text(x=c(8), y=rep(13, 1), label="*", cex=1.5)

segments(x0=6,x1=8, y0=10, y1=10, lwd=2, col="gray")
segments(x0=6,x1=6, y0=10, y1=9, lwd=2)
text(x=c(7), y=rep(10, 1), label="*", cex=1.5)

segments(x0=7,x1=8, y0=7, y1=7, lwd=2, col="gray")
segments(x0=7,x1=7, y0=7, y1=6, lwd=2)

plot.new()
par(lheight=0.7)
text(x=0.3, y=0.5, label="Number of tumors", cex=1.0, srt=90, xpd=TRUE)
text(label="C", x=-0.8, y=1.03, font=2, cex=2, xpd=TRUE)



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
        pch=pch[ordered_treatments[i]], jitter=0.25, cex=1.0, add=T)
    }

    text(x=0.35, y=1.15*max(rabunds), label=otu_label, pos=4, font=2, cex=0.8, xpd=TRUE)
}


otu_plot(baseline_rabund[,ordered_otus[1]], ordered_treatments, tax[ordered_otus[1],"Taxonomy"])
otu_plot(baseline_rabund[,ordered_otus[2]], ordered_treatments, tax[ordered_otus[2],"Taxonomy"])
otu_plot(baseline_rabund[,ordered_otus[3]], ordered_treatments, tax[ordered_otus[3],"Taxonomy"])
otu_plot(baseline_rabund[,ordered_otus[4]], ordered_treatments, tax[ordered_otus[4],"Taxonomy"])


plot.new()
text(x=0.3, y=0.5, label="Relative abundance (%)", cex=1.0, srt=90, xpd=TRUE)
text(label="D", x=-0.8, y=1.03, font=2, cex=2, xpd=TRUE)




plot.new()
text(x=seq(0.05,1,length=8), y=1.2, labels[ordered_treatments], srt=70, cex=0.8, font=1, pos=2, xpd=TRUE)



plot.new()
par(mar=c(0,0,0,0))
img <- readPNG("results/pictures/tumor_images.png")
rasterImage(img, -0.1, 0.0, 1.0, 0.95)
text(label="B", x=-0.1, y=1.03, font=2.5, cex=2, xpd=TRUE)



nmds_file <- "data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.nmds.axes"
nmds <- read.table(file=nmds_file, header=T, stringsAsFactors=FALSE)
treatment_group <- gsub("^([^_]*)_.*", "\\1", nmds$group)
par(mar=c(4,4,0.5,0.5))
plot(nmds$axis2~nmds$axis1, xlab="NMDS Axis 1", ylab="NMDS Axis 2", col=clrs[treatment_group],pch=pch[treatment_group], ylim=c(-0.8,0.8), axes=FALSE)
box()
axis(1, label=format(seq(-0.6,0.6,0.3), nsmall=1), at=seq(-0.6,0.6,0.3))
axis(2, label=format(seq(-0.6,0.6,0.3), nsmall=1), at=seq(-0.6,0.6,0.3), las=1)

text(label="E", x=-0.88, y=0.83, font=2.5, cex=2, xpd=TRUE)

dev.off()
