library("png")
source("code/rf_baseline_analysis.R")

design <- c(1,1,10,
            6,2,10,
            6,3,10,
            6,4,10,
            6,5,10,
            8,7,10,
            0,9,10)

design <- matrix(design, ncol=3, byrow=TRUE)

z <- layout(design, width=c(0.2, 1, 1), height=c(2,1,1,1,1,2,1.2))


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
    Metro = paste0("Metronidozole (N=",counts["Metro"] , ")"),
    VancStrep = paste0("\U0394Metronidozole (N=", counts["Metro"], ")"),
    Strep = paste0("Streptomycin (N=",counts["Strep"] , ")"),
    VancMetro = paste0("\U0394Streptomycin (N=",counts["Strep"] , ")"),
    Vanc = paste0("Vancomycin", " (N=",counts["Vanc"] , ")"),
    StrepMetro = paste0("\U0394Vancomycin (N=",counts["StrepMetro"], ")") )


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
text(label="B", x=-0.8, y=1.03, font=2, cex=2, xpd=TRUE)


plot(NA, xlim=c(0.5,8.5), xlab="", ylim=c(0,25), axes=F, cex.lab=1.4)
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

plot.new()
par(lheight=0.7)
text(x=0.2, y=0.5, label="Number of tumors", cex=1.0, srt=90, xpd=TRUE)
text(label="C", x=-0.8, y=1.03, font=2, cex=2, xpd=TRUE)


plot.new()
text(x=seq(0.05,1,length=8), y=1.2, labels[ordered_treatments], srt=70, cex=0.8, font=1, pos=2, xpd=TRUE)



plot.new()
img <- readPNG("results/pictures/tumor_images.png")
lim <- par()
rasterImage(img, -0.1, 0.5, 1.0, 1.0)
text(label="D", x=0.0, y=1.03, font=2.5, cex=2, xpd=TRUE)
