#read in a mothur shared file that has already been subsampled

read_shared_file <- function(shared_file, samples_to_remove){

    path_to_file_name <- paste0(shared_file)
    #the second column has the group names and there is a header
    shared <- read.table(file=path_to_file_name, row.names=2, header=T, stringsAsFactors=FALSE)

    #get rid of the first two columns, which are not needed - these were
    #the first and third columns in the original shared file
    cols_to_remove <- c("label", "numOtus")
    shared <- shared[, !colnames(shared)%in% cols_to_remove]

    #remove those rows that are labeled as being samples that only received DSS
    rows_to_remove <- grepl(samples_to_remove, rownames(shared))
    shared <- shared[!rows_to_remove, ]
}


# given a table of counts where the rows are the samples and the columns are the
# otus, output the relative abundance for each OTU in each sample

get_rel_abund <- function(file_name, samples_to_remove){

    shared_table <- read_shared_file(file_name, samples_to_remove)
    #calculate the number of sequences in each sample: assumes every sample
    #has the same number of sequences
    n_seqs <- sum(shared_table[1,])

    #calculate the relative abundances
    shared_table/n_seqs
}



#The sample names go treatment_eartag_day, so let's split them on the '_'
#and then return the appropriate slot in the resulting vector
get_name_part <- function(sample_name, column){
    name_parts <- unlist(strsplit(sample_name, "_"))
    return(name_parts[column])
}

#treatment is in slot 1
get_treatment <- function(sample_name){
    get_name_part(sample_name, column=1)
}

#eartag is in slot 2
get_eartag <- function(sample_name){
    get_name_part(sample_name, column=2)
}

#eartag is in slot 3; neg indicates that the sample was prior to day 0;
#need to convert string to a number
get_day <- function(sample_name){
    text.day <- get_name_part(sample_name, column=3)
    text.day <- gsub("d", "", text.day)
    text.day <- gsub("neg", "-", text.day)
    as.numeric(text.day)
}



#keep those OTUs that are more than threshold across all samples
filter_rabund <- function(rabund, threshold){
    filter <- apply(rabund, 2, mean) > threshold
    rabund[,filter]
}


#build random forest with the filtered rabund table; dependent would be
#tumor_counts
get_forest <- function(rabund, dependent, threshold, n_trees=10000){

    #keep those OTUs that are more than threshold across all samples
    rabund_filtered <- filter_rabund(rabund, threshold)

    #build random forest model where we predict `dependent` based on the
    #filtered rabund data using `n_trees`
    randomForest(dependent ~ ., data=rabund_filtered, importance=TRUE, ntree=n_trees)
}


#extract the rsquared value from the forest
get_rsq <- function(forest, n_trees=10000){
    forest$rsq[n_trees]
}


#extract the number of features used to build the forest
get_n_otus <- function(forest){
    nrow(forest$importance)
}


#see what the rsquared value looks like for forests that are built stepwise
simplify_model <- function(dependent, forest, rabund){
    #gini coefficient is in the second column of the importance data frame
    importance <- importance(forest)
    sorted_importance <- importance[order(importance[,2], decreasing=T),]

    notus <- nrow(importance)
    rf_simplify_rsq <- rep(0, notus)

    #add each successive OTU's data and rebuild the model; save the Rsquared value
    #can't have a model with only one feature; start at i = 2.
    for(i in 2:notus){

        #extract the appropriate columns
        simple_rabund <- rabund[,colnames(rabund) %in% rownames(sorted_importance)[1:i]]
        #build the model
        rf_simplify <- randomForest(dependent ~ ., data=simple_rabund,
                                                importance=TRUE, ntree=n_trees)
        #get the Rsquared values
        rf_simplify_rsq[i] <- rf_simplify$rsq[n_trees] #percent variance explained
    }
    return(rf_simplify_rsq)
}


#information for plotting...
clrs <- c(Metro = "red",
    VancStrep = "red",
    Strep = "darkgreen",
    VancMetro = "darkgreen",
    Vanc = "blue",
    StrepMetro = "blue",
    AllAbs = "black",
    NoAbs = "black")

#single Abx get a solid dot, the deletion gets an open dot
pch <- c(Metro = 19,
    VancStrep = 1,
    Strep = 19,
    VancMetro = 1,
    Vanc = 19,
    StrepMetro = 1,
    AllAbs = 19,
    NoAbs = 1)

#these are needed for the legend
labels <- c(NoAbs = "No antibiotics",
    AllAbs = "All antibiotics",
    Metro = "Metronidozole",
    VancStrep = "\u0394 Metronidozole",
    Strep = "Streptomycin",
    VancMetro = "\u0394 Streptomycin",
    Vanc = "Vancomycin",
    StrepMetro = "\u0394 Vancomycin" )


#forest <- rf_baseline_forest
plot_importance <- function(forest){

    par(mar=c(3,10,0.5,0.5))
    #read in the taxonomy file
    tax <- read.table(file="data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy", header=T, row.names=1)

    #the next three lines extract the last named taxonomic level for each OTU
    tax$Taxonomy <- gsub("unclassified.*", "", tax$Taxonomy)
    tax$Taxonomy <- gsub("\\(\\d*\\);$", "", tax$Taxonomy)
    tax$Taxonomy <- gsub(".*;", "", tax$Taxonomy)
    tax$Taxonomy <- gsub("_", " ", tax$Taxonomy)
    tax$Taxonomy <- gsub(" sensu stricto", " ", tax$Taxonomy)
    otu_tax_labels <- paste0(tax[,2]," (",gsub("tu0*", "TU ", rownames(tax)) ,")")
    names(otu_tax_labels) <- rownames(tax)

    importance <- importance(forest)
    sorted_importance <- importance[order(importance[,2], decreasing=T),2]

    plot(NA, yaxt="n", xlab="", ylab="",
        xlim=c(0, max(sorted_importance)),
        ylim=c(1, length(sorted_importance)))
    abline(h=1:length(sorted_importance), lty=3, col="gray")
    points(x=rev(sorted_importance), y=1:length(sorted_importance), pch=19)
    mtext(side=2, line=8.5, adj=0, at=1:length(sorted_importance), text=rev(otu_tax_labels[names(sorted_importance)]), las=2, cex=0.7)
    mtext(side=1, text="Gini coefficient", line=2.0)

}


#fit the full model back to the original data
plot_forest_fit <- function(observed, forest, rabund, treatment){

    par(mar=c(3,3,0.5,0.5))
    forest_fit <- predict(forest, rabund)

    #want the x and y-axes to have the same limits, so let's find the largest
    #tumor value among the observed and predicted
    max_value <- max(c(forest_fit, observed))

    labels_no_n <- gsub(" \\(.*\\)", "", labels)

    plot(forest_fit~observed, pch=pch[treatment],
        col=clrs[treatment], ylim=c(0,max_value), xlim=c(0,max_value),
        cex=1, xlab="", ylab="", yaxt="n")
    legend(x=-1, y=max_value*1.05, legend=labels_no_n, pch=pch[names(labels)],
            col=clrs[names(labels)], cex=0.8, bty="n")
    mtext(side=1, text="Observed number of tumors", line=2.0)
    mtext(side=2, text="Predicted number of tumors", line=2.0)
    axis(2, las=2)

}

#Plot top features' relative abundance versus the tumor counts for the mice that
#they came from...
plot_top_features <- function(tumor_counts, sorted_importance, rabund, treatment){
    #read in the taxonomy file
    tax <- read.table(file="data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy", header=T, row.names=1)

    #the next three lines extract the last named taxonomic level for each OTU
    tax$Taxonomy <- gsub("unclassified.*", "", tax$Taxonomy)
    tax$Taxonomy <- gsub("\\(\\d*\\);$", "", tax$Taxonomy)
    tax$Taxonomy <- gsub(".*;", "", tax$Taxonomy)

    #let's just use the top six OTUs, based on our inspection of the importance plot
    #and combine the OTU name with its taxonomy and Gini coefficient
    otus <- rownames(sorted_importance)[1:6]
    pretty_otus <- gsub("Otu0*", "OTU ", otus)
    otu_labels <- paste0("(", pretty_otus, ")", "\nGini: ", round(sorted_importance[1:6,2], 1))
    otu_labels <- paste(tax[otus,2], otu_labels, sep=" ")


    par(mar=c(0.5,0.5,0.5,0.5))

    design <- matrix(1:6, nrow=2, byrow=T)
    design <- cbind(c(7,7), design)
    design <- rbind(design, c(0,8,8,8))
    layout(design, widths=c(0.3,1,1,1), heights=c(1,1,0.3))

    for(i in 1:6){

        #get the row and column number for each spot in the layout
        row <- ifelse(i<=3,1,2)
        column <- ifelse(i<=3, i, i-3)

        #extract the relative abundance data for this OTU
        rel_abund <- rabund[,otus[i]]

        #plot the relative abundance with the number of tumors for each animal. plot
        #on consistent log scaled x-axis for all OTUs. will throw errors because it
        #can't plot zeroes on a log scale
        plot(rel_abund,tumor_counts, log="x", pch=pch[treatment],
            col=clrs[treatment], ylab="", xlab="",
            xlim=c(1e-4, 1), ylim=c(0,27), yaxt="n",
            xaxt="n", cex.lab=1.5)

        #want to plot the number of tumors for those mice that had a zero relative
        #abundance
        zeroes <- rel_abund == 0

        #jitter the points so that they don't fall on top of each other
        x_zeroes <- runif(sum(zeroes),1.0e-4,1.5e-4)#rep(1.2e-4, sum(zeroes))
        points(x=x_zeroes, tumor_counts[zeroes], pch=pch[treatment[zeroes]],
                col=clrs[treatment[zeroes]])

        #create a vertical line to denote the limit of detection
        abline(v=2.2e-4, col="gray")

        #put the OTU label in the upper left corner of the plot
        text(x=0.8e-4, y=25, label=otu_labels[i], pos=4, font=2)

        #if it's on the bottom row, put a customized axis indicating the % rabund
        if(row == 2){
            axis(1, at=c(1.25e-4, 1e-3,1e-2,1e-1,1),
                    label=c("0", "0.1", "1", "10", "100"),
                    cex.axis=1.5)
        }

        #if it's in the first column turn the axis labels to be horizontal
        if(column == 1){
            axis(2, las=2, cex.axis=1.5)
        }
    }

    plot.new()
    text(x=0.15, y=0.5, label="Number of tumors", cex=2, srt=90)

    plot.new()
    text(x=0.5, y=0.2, label="Relative abundance (%)", cex=2)
}
