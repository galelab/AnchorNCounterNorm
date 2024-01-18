

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
#' @import reshape2 norm_matrix_melt <- melt(norm_matrix)
#' @examples
#' data <- hk_gene_stats(counts)


ratio_normalization <- function (data, hkgenes, output_dir=getwd(), save.fig=TRUE, separte_norm_counts_by_condition=NULL) {
    countsall <- data$counts[data$counts$Class%in% c("Endogenous", "Housekeeping"), 3:length(colnames(data$counts))]
    rownames(countsall) <- data$counts[data$counts$Class%in% c("Endogenous", "Housekeeping"), "Gene_Name"]
    hkcountgeomeans <- apply(countsall[hkgenes, ], MARGIN = 2, geoMean)
    print(hkcountgeomeans)
    norm_ratio <- data$counts[data$counts$Class%in% c("Endogenous"), 3:length(colnames(data$counts))]
    rownames(norm_ratio) <-  data$counts[data$counts$Class%in% c("Endogenous"),  "Gene_Name"]
    for (col in colnames(norm_ratio)) {
        norm_ratio[, col] <- norm_ratio[, col] / hkcountgeomeans[col]
    }
    pl <- boxplot_normalization(norm_ratio, data$meta.data, output_dir, save.fig = save.fig)
    print(head(norm_ratio))
    write.csv(log2(norm_ratio), file.path(output_dir, "lognormratiocounts.csv"))
    return(list("counts_ratio"=norm_ratio, "log_counts_ratio"=log(norm_ratio), "boxplot"=pl))
}

boxplot_normalization <- function(norm_matrix ,meta.data, output_dir=getwd(), save.fig=TRUE, separte_norm_counts_by_condition=NULL) {
    norm_matrix_melt <- melt(norm_matrix)
    if (length(unique(norm_matrix_melt$variable)) < 20) {
        fontsize <- 6
    } else {
        fontsize <- 3
    }
    pl <- ggplot(norm_matrix_melt, aes(x = variable, y=log2(value))) +
        geom_boxplot(fill = "#ff615d") + theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = fontsize), 
        axis.title.x=element_blank())
    if(isTRUE(save.fig)) {
        message("STATUS: saving figure as png file here: ", output_dir)
        ggsave(file.path(output_dir, "boxplotnormlog2.png"), width = 4.5, height = 3, bg = "white", dpi = 300)
        ggsave(file.path(output_dir, "boxplotnormlog2.pdf"), width = 4.5, height = 3, bg = "white", dpi = 300)

    }
    return(pl)
}