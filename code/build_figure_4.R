source("code/rf_baseline_analysis.R")


get_dist_lt <- function(file_name){

    dist_vector <- scan(file_name, what="", sep="\n")

    n_samples <- as.numeric(dist_vector[1])
    dist_vector <- dist_vector[-1]

    sample_names <- character(n_samples)
    distance_matrix <- matrix(rep(0, n_samples*n_samples), nrow=n_samples)

    for(i in 1:n_samples){

        split_line <- unlist(strsplit(dist_vector[i], "\t"))
        sample_names[i] <- split_line[1]

        if(i > 1){
            distance_matrix[i,1:(i-1)] <- as.numeric(split_line[2:i])
        }
    }

    distance_matrix <- distance_matrix + t(distance_matrix)
    rownames(distance_matrix) <- sample_names
    colnames(distance_matrix) <- sample_names
    distance_matrix
}


get_start_end_dist <- function(u_treat_et, days){
    subset <- days[treat_et == u_treat_et]
    start <- names(subset[subset==0])
    end <- names(sort(subset, decreasing=TRUE))[1]
    distances[start, end]
}


dist_file <- "data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.dist"

distances <- get_dist_lt(dist_file)
treat_et <- gsub("_d.*", "", rownames(distances))
days <- sapply(rownames(distances), get_day)
start_end_dist <- unlist(sapply(unique(treat_et), get_start_end_dist, days))

tumor_count.file <- "data/process/tumor_counts.tsv"
tumor_counts <- read.table(file=tumor_count.file, header=T)
tumor_counts_prepost <- tumor_counts[names(start_end_dist), "n_tumors"]
dist_tumor_cor <- cor.test(tumor_counts_prepost, start_end_dist,
                                method="spearman")

save(dist_tumor_cor, file="data/process/distance_tumor_correlation.Rdata")

pdf(file="results/figures/figure_4.pdf", width=3.75, height=7)
design <- matrix(c(2,1,
                    5,3,
                    5,4,
                    0,6), byrow=T, nrow=4)
z <- layout(design, width=c(0.15,1), height=c(1.5,1,1,0.3))

par(mar=c(3,0.5,0.5,0.5))
treatments <- gsub("_.*", "", names(start_end_dist))
plot(tumor_counts_prepost~start_end_dist, ylim=c(0,25), xlim=c(0,1),
    col=clrs[treatments],
    pch=pch[treatments], axes=FALSE, xlab="", ylab="")
mtext(side=1, text="Distance between communities on days 0 and 74", line=2, cex=0.7)
axis(1)
axis(2, las=2)
box()

#labels_no_n <- gsub(" \\(.*\\)", "", labels)
#legend(x=-0.025, y=26, legend=labels_no_n, pch=pch[names(labels)],
#        col=clrs[names(labels)], cex=0.9, bty="n")

par(mar=c(0,0,0,0))
plot(NA, xlim=c(0,1), ylim=c(0,1), axes=F, xlab="", ylab="", xpd=TRUE)
text(x=0.3, y=0.58, label="Number of tumors", srt=90, cex=1.5)
text(x=0.3, y=0.98, label="A", srt=0, cex=2, font=2)


par(mar=c(0.5,0.5,0.5,0.5))
source("code/make_timeline_plot.R")

par(mar=c(0,0,0,0))
plot(NA, xlim=c(0,1), ylim=c(0,1), axes=F, xlab="", ylab="", xpd=TRUE)
text(x=0.3, y=0.5, label="Relative abundance (%)", srt=90, cex=1.5)
text(x=0.3, y=1, label="B", srt=0, cex=2, font=2)

plot(NA, xlim=c(0.1,1), ylim=c(0,1), axes=F, xlab="", ylab="")
text(x=0.55, y=0.25, label="Days since AOM administration")
dev.off()
