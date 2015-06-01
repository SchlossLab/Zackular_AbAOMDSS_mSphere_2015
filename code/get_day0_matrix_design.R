###############################################################################
#
#   get_day0_matrix.R
#
#   Create a sub distance matrix and design file for day 0 from the full
#   distance matrix. Assumes that the sample names contain the day number at
#   the end of the sample name using a "_dXX" naming convention. Outputs a
#   square distance matrix
#
###############################################################################

get_dist_lt <- function(file_name){

    dist_vector <- scan(file_name, what="", sep="\n", quiet=TRUE)

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

get_day_samples <- function(dist_lt, day){

    sample_ids <- rownames(dist_lt)
    right_day <- grepl(paste0("_d", day, "$"), sample_ids)
    dist_lt[right_day, right_day]

}

get_design <- function(dist_lt){

    sample_ids <- rownames(dist_lt)
    groups <- gsub("^([^_]*)_.*", "\\1", sample_ids)

    cbind(sample_ids, groups)

}


phylip_filename <- "data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.dist"
sub_dist_filename <- "data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.dist"
design_filename <- "data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.design"


dist_lt <- get_dist_lt(phylip_filename)
day0_dist <- get_day_samples(dist_lt, 0)

write(file=sub_dist_filename, nrow(day0_dist))
write.table(file=sub_dist_filename, day0_dist, quote=F, col.names=F, sep="\t", append=TRUE)


day0_design <- get_design(day0_dist)
write.table(file=design_filename, day0_design, quote=F, col.names=F, row.names=F, sep="\t")
