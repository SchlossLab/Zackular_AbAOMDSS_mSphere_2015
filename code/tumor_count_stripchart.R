#build plot that shows the number of tumors that the mice received who were
#given the various antibiotic treatments.

#get some of the functions/colors/etc. from the rf file
source("code/rf_baseline_analysis.R")

tumors <- read.table(file="data/process/tumor_counts.tsv", header=T, row.names=1)
treatments <- sapply(rownames(tumors), get_treatment)

#get the mean number of tumors in each treatment
summarize_counts <- aggregate(tumors$n_tumors, by=list(treatments), mean)
means <- summarize_counts$x
names(means) <- summarize_counts$Group.1


#order the treatments by the number of tumors
decreasing_order <- order(means, decreasing=T)
ordered_treatments <- levels(as.factor(treatments))[decreasing_order]

counts <- table(treatments)

labels <- c(NoAbs = paste0("No antibiotics (N=",counts["NoAbs"],")"),
    AllAbs = paste0("All antibiotics (N=",counts["AllAbs"] , ")"),
    Metro = paste0("Metronidozole (N=",counts["Metro"] , ")"),
    VancStrep = paste0("\U0394Metronidozole (N=", counts["Metro"], ")"),
    Strep = paste0("Streptomycin (N=",counts["Strep"] , ")"),
    VancMetro = paste0("\U0394Streptomycin (N=",counts["Strep"] , ")"),
    Vanc = paste0("Vancomycin", " (N=",counts["Vanc"] , ")"),
    StrepMetro = paste0("\U0394Vancomycin (N=",counts["StrepMetro"], ")") )


cairo_pdf(file="tumor_counds.pdf", width=4.5,5.5)
par(mar=c(10, 5, 0.5, 0.5))

plot(1, xlim=c(0.5,8.5), xlab="", ylim=c(0,max(tumors$n_tumors)), axes=F,
    ylab="Number of tumors", cex.lab=1.4)
segments(x0=seq(1:8)-0.4, x1=seq(1:8)+0.4, y0=means[decreasing_order], y1=means[decreasing_order], lwd=3)
box()
axis(2, las=2)
axis(1, at=1:8, label=labels[ordered_treatments], las=2, tick=F)


#build a stripchart
for(i in 1:length(ordered_treatments)){
    stripchart(at=i, x=tumors[treatments==ordered_treatments[i], "n_tumors"],
    vertical=T, method="jitter", col=clrs[ordered_treatments[i]],
    pch=pch[ordered_treatments[i]], jitter=0.25, cex=1.5, add=T)
}


dev.off()
