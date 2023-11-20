

#' ratio normalization
#'
#' This implements ratio normalization
#' @param counts  matrix counts
#' @param hkgenes list of hk genes to use to get ratio
#' @param output_dir directory for results (default is working directory)
#' @keywords house keeping genes
#' @export
#' @import EnvStats
#' @import ggplot2
#' @examples
#' data <- hk_gene_stats(counts)


ratio_normalization <- function (data, hkgenes, output_dir=getwd()) {
    countsall <- data$counts[data$counts$Class%in% c("Endogenous", "Housekeeping"), 3:length(colnames(data$counts))]
    rownames(countsall) <- data$counts[data$counts$Class%in% c("Endogenous", "Housekeeping"), "Gene_Name"]
    hkcountgeomeans <- apply(countsall[hkgenes, ], MARGIN = 2, geoMean)
    print(hkcountgeomeans)
    norm_ratio <- data$counts[data$counts$Class%in% c("Endogenous"), 3:length(colnames(data$counts))]
    rownames(norm_ratio) <-  data$counts[data$counts$Class%in% c("Endogenous"),  "Gene_Name"]
    for (col in colnames(norm_ratio)) {
        norm_ratio[, col] <- norm_ratio[, col] / hkcountgeomeans[col]
    }
    pl <- boxplot_normalization(norm_ratio, output_dir)
    print(head(norm_ratio))
    write.csv(log2(norm_ratio), file.path(output_dir, "lognormratiocounts.csv"))
    return(list("counts_ratio"=norm_ratio, "log_counts_ratio"=log(norm_ratio), "boxplot"=pl))
}

boxplot_normalization <- function(norm_matrix, output_dir=getwd()) {
    png(file.path(output_dir, "boxplotnormlog2.png"), width=6, height=5, units="in", res=100)
    pl <-boxplot(log2(norm_matrix),
        col = "#ff615d", main = "log norm ratio counts", cex.axis = .6, las = 2,
        xlab = "",frame = FALSE, ylab = "log(Norm Ratio Counts)"
    )
    dev.off()
    return(pl)
}