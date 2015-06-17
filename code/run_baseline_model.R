source("code/rf_baseline_analysis.R")
library("randomForest")

set.seed(19760620)

shared_file <- "data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared"
rel_abund <- get_rel_abund(shared_file, "DSS")

treatment <- sapply(rownames(rel_abund), get_treatment)
eartag <- sapply(rownames(rel_abund), get_eartag)
day <- sapply(rownames(rel_abund), get_day)

baseline <- day == 0                        #baseline samples are from day 0
treatment_baseline <- treatment[baseline]   #update the metadata
eartag_baseline <- eartag[baseline]         #update the metadata
day_baseline <- day[baseline]               #update the metadata
rabund_baseline <- rel_abund[baseline,]     #extract the baseline samples

#what mice do we have baseline data for? necessary because the tumor data
#is indexed by `treatment_eartag`
mice_baseline <- paste(treatment_baseline, eartag_baseline, sep="_")

tumor_count.file <- "data/process/tumor_counts.tsv"
tumor_counts <- read.table(file=tumor_count.file, header=T)
tumor_counts_baseline <- tumor_counts[mice_baseline,]   #get the tumor count data


n_trees <- 10000 #number of trees for random forest

rabund_baseline_filtered <- rabund_baseline[,apply(rabund_baseline, 2, sum)>0]

n_full_features <- ncol(rabund_baseline_filtered)


rf_baseline_forest <- get_forest(rabund=rabund_baseline_filtered,
							dependent=tumor_counts_baseline, n_trees=n_trees)

rf_simplify_rsq <- simplify_model(tumor_counts_baseline, rf_baseline_forest,
    						rabund_baseline_filtered, 15)

n_features <- which.max(rf_simplify_rsq)

decrease_mse <- importance(rf_baseline_forest,1)
feature_order <- order(decrease_mse, decreasing=TRUE)
top_features <- names(decrease_mse[feature_order,])[1:n_features]

rabund_baseline_top_features <- rabund_baseline[,as.character(top_features)]
rf_baseline_top_features_forest <- get_forest(rabund=rabund_baseline_top_features,
							dependent=tumor_counts_baseline, n_trees=n_trees)


tax <- read.table(file="data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy", header=T, row.names=1)
tax$Taxonomy <- gsub("unclassified.*", "", tax$Taxonomy)
tax$Taxonomy <- gsub("\\(\\d*\\);$", "", tax$Taxonomy)
tax$Taxonomy <- gsub(".*;", "", tax$Taxonomy)

baseline_all_rsq <- rf_baseline_forest$rsq[n_trees]
baseline_top_rsq <- rf_baseline_top_features_forest$rsq[n_trees]
baseline_features <- tax[top_features,]

save(n_full_features, baseline_all_rsq, baseline_top_rsq, baseline_features, file="data/process/baseline_model.Rdata")


cairo_pdf(file="results/figures/figure_S1.pdf", width=5.0, height=5.0)
	plot(rf_simplify_rsq, xlim=c(2,15), ylim=c(0.5, 0.7), xlab="Number of features",
		ylab="Variance explained (%)", pch=19)
	text(x=n_features, y=0.70, "*", cex=2)
dev.off()



cairo_pdf(file="results/figures/figure_2.pdf", width=3.5, height=3.5)
	plot_forest_fit(tumor_counts_baseline, rf_baseline_top_features_forest,
					rabund_baseline_top_features, treatment_baseline)
dev.off()



cairo_pdf(file="results/figures/figure_3.pdf", width=7.5, height=10)
	plot_baseline_features(tumor_counts_baseline, rf_baseline_top_features_forest,
						rabund_baseline_top_features, treatment_baseline)
dev.off()
