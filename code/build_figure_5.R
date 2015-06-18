library("png")

pdf("results/figures/figure_5.pdf", width=3.5, height=3.5)

abx <- read.table(file="data/process/tumor_counts.tsv", header=T)
intervention <- read.table(file="data/process/intervention_tumor_counts.tsv", header=T)

layout(matrix(c(1,2), byrow=T, nrow=1), width=c(1.5,1.5), height=1)

par(mar=c(6,4,0.5,0.25))

plot(NA, xlim=c(0.5,4.5), ylim=c(0,25), xlab="", ylab="Number of tumors", axes=F)
stripchart(at=1, x=abx[grepl("NoAbs", rownames(abx)), "n_tumors"],
    vertical=T, method="jitter", col="black", pch=1, jitter=0.15, cex=1.0, add=T)
stripchart(at=2, x=abx[grepl("AllAbs", rownames(abx)), "n_tumors"],
    vertical=T, method="jitter", col="black", pch=19, jitter=0.15, cex=1.0, add=T)
stripchart(at=3, x=intervention[intervention$experiment=="Intervention1", "n_tumors"],
    vertical=T, method="jitter", col="black", pch=18, jitter=0.15, cex=1.5, add=T)
stripchart(at=4, x=intervention[intervention$experiment=="Intervention2", "n_tumors"],
    vertical=T, method="jitter", col="black", pch=17, jitter=0.15, cex=1.5, add=T)

text(y=-9, x=(1:4)-0.6, label=c("No antibiotics", "All antibiotics", "Intervention 2", "Intervention 3"), srt=70, cex=0.8, xpd=TRUE)
axis(2, las=1)
axis(1, tick=TRUE, labels=FALSE)
box()
mtext(text="A", side=2, line=2.5, at=25, font=2.5, cex=2, las=2, xpd=TRUE)


par(mar=c(6,1.5,0.0,0.0))
plot(c(0, 150), c(0, 150), type = "n", xlab = "", ylab = "", axes=FALSE)
img <- readPNG("results/pictures/tumor_images_intervention.png")
lim <- par()
rasterImage(img, -10, -10, 160, 160, xpd=TRUE)
text(label="B", x=-30.0, y=143, font=2.5, cex=2, xpd=TRUE)

dev.off()
