source("code/rf_baseline_analysis.R")
library("randomForest")

set.seed(19760620)

shared_file <- "data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared"
rel_abund <- get_rel_abund(shared_file, "DSS")

treatment <- sapply(rownames(rel_abund), get_treatment)
eartag <- sapply(rownames(rel_abund), get_eartag)
day <- sapply(rownames(rel_abund), get_day)


get_end_sample <- function(treat_et){
    treat_et_days <- day[grep(treat_et, names(day))]
    names(which.max(treat_et_days))
}

unique_treat_et <- unique(paste(treatment, eartag, sep="_"))
ends <- sapply(unique_treat_et, get_end_sample)

treatment_final <- treatment[ends]   #update the metadata
eartag_final <- eartag[ends]         #update the metadata
day_ends <- day[ends]               #update the metadata
rabund_final <- rel_abund[ends, ]

tumor_count.file <- "data/process/tumor_counts.tsv"
tumor_counts <- read.table(file=tumor_count.file, header=T)
tumor_counts_final <- tumor_counts[unique_treat_et, ]

n_trees <- 10000 #number of trees for random forest
rabund_final_filtered <- rabund_final[,apply(rabund_final, 2, sum)>0]

rf_final_forest <- get_forest(rabund=rabund_final_filtered,
							dependent=tumor_counts_final, n_trees=n_trees)

rf_simplify_rsq <- simplify_model(tumor_counts_final, rf_final_forest,
    						rabund_final_filtered, 15)

n_features <- which.max(rf_simplify_rsq)



decrease_mse <- importance(rf_final_forest,1)
feature_order <- order(decrease_mse, decreasing=TRUE)
top_features <- names(decrease_mse[feature_order,])[1:n_features]

rabund_final_top_features <- rabund_final[,as.character(top_features)]
rf_final_top_features_forest <- get_forest(rabund=rabund_final_top_features,
							dependent=tumor_counts_final, n_trees=n_trees)


tax <- read.table(file="data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy", header=T, row.names=1)
tax$Taxonomy <- gsub("unclassified.*", "", tax$Taxonomy)
tax$Taxonomy <- gsub("\\(\\d*\\);$", "", tax$Taxonomy)
tax$Taxonomy <- gsub(".*;", "", tax$Taxonomy)

final_all_rsq <- rf_final_forest$rsq[n_trees]
final_top_rsq <- rf_final_top_features_forest$rsq[n_trees]
final_features <- tax[top_features,]


save(final_all_rsq, final_top_rsq, final_features, file="data/process/final_model.Rdata")

cairo_pdf(file="results/figures/figure_S2.pdf", width=5.0, height=5.0)
	plot(rf_simplify_rsq, xlim=c(2,15), ylim=c(0.5, 0.7), xlab="Number of features",
		ylab="Variance explained (%)", pch=19)
	text(x=n_features, y=0.70, "*", cex=2)
dev.off()

cairo_pdf(file="results/figures/figure_S3.pdf", width=3.5, height=3.5)
	plot_forest_fit(tumor_counts_final, rf_final_top_features_forest,
					rabund_final_top_features, treatment_final)
dev.off()


cairo_pdf(file="results/figures/figure_S4.pdf", width=7.5, height=10)
	plot_final_features(tumor_counts_final, rf_final_top_features_forest,
						rabund_final_top_features, treatment_final)
dev.off()
