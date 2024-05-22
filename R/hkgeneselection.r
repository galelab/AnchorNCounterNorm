
#' house keeping gene analysis
#'
#' This assesses which genes should be used for ratio normalization 
#' @param data  counts compiled from rcc files in load_nCounter_files step 
#' @param group.by what group counts of each house keeping gene, individual samples are recomme5ded
#' @param min_num_hk_genes minimum number of house keeping genes to keep for analysis
#' @param manually_selected_hk_genes list of house keeping genes manually selected 
#' @param output_dir directory for results (default is working directory)
#' @keywords house keeping genes
#' @export
#' @import EnvStats
#' @import ggplot2
#' @import reshape2
#' @import ggthemes
#' @examples
#' data <- hk_gene_stats(counts)
#' 
#'

hk_gene_stats <- function(data, group.by="sample",
    min_num_hk_genes = 5, manually_selected_hk_genes = FALSE, output_dir = getwd()) {
    get_CV_4_alltypes(data$counts, output_dir = output_dir)
    if (isFALSE(manually_selected_hk_genes)) {
        lowestcv <- identifyhkgenes(data$counts, hknum = min_num_hk_genes, output_dir = output_dir)
        pl <- violin_plots_hkgenes(data$counts[data$counts$Gene_Name %in% lowestcv$hkgenes, ], 
            data$meta.data, group.by=group.by,
            output_dir = output_dir, "lowestCVgenesviolinplot")
        pl2 <- violin_plots_hkgenes(data$counts[data$counts$Class =="Housekeeping", ], 
            data$meta.data, group.by=group.by,
            output_dir = output_dir, "HKgenesCVviolinplot")
        hkstats <- get_CV_4_HKgenes(data$counts[data$counts$Class =="Housekeeping", ], output_dir = output_dir)
        return(list("lowestCVstats" = lowestcv, "hkstats" = hkstats$hkstats, "lowestCVviolinplot" = pl, "HKgenesCVviolinplot" = pl2))

    } else {
        pl <- violin_plots_hkgenes(data$counts[data$counts$Gene_Name %in% manually_selected_hk_genes, ], 
            data$meta.data, group.by=group.by,
            output_dir = output_dir, "hkviolinplotmanually")
        return(list("violinplot" = pl))
    }
}

calc_CV <- function(x) {
    sd(x) / mean(x)
}

get_CV_4_alltypes <- function(df, output_dir) {
    dfdata <- df[, 3:length(colnames(df))]
    gene_names <- df[, "Gene_Name"]
    class <- df[, "Class"]
    print(length(class))
    dfdata <- as.matrix(as.data.frame(dfdata))
    class(dfdata) <- "numeric"
    print(dim(dfdata))
    CV_all <- apply(dfdata, MARGIN = 1, calc_CV)
    mean_all <- apply(dfdata, MARGIN = 1, mean)
    sd_all <- apply(dfdata, MARGIN = 1, sd)
    geomean_all <- apply(dfdata, MARGIN = 1, geoMean)
    message("STATUS: Loading data into data.frame")
    alldf <- data.frame(matrix(nrow = length(rownames(dfdata)), ncol = 6))
    colnames(alldf) <- c("mean", "geomean", "sd", "CV", "genes", "Class")
    alldf$mean <- as.numeric(mean_all)
    alldf$geomean <- as.numeric(geomean_all)
    alldf$sd <- as.numeric(sd_all)
    alldf$CV <- as.numeric(CV_all)
    alldf$Class <- class
    alldf$genes <- gene_names
    alldf <- alldf[order(CV_all), ]
    write.csv(alldf, file.path(output_dir, "allCVagenesstats.csv"), quote = F, row.names = F)
    ggplot(alldf, aes(x=Class, y=CV_all, color=Class)) + geom_boxplot() +theme_Publication(base_size=8) + 
        scale_fill_manual(values = c("Endogenous"="#ff615d", "Housekeeping"="#004c7a", "Negative"="#ffd400", "Positive"="#61c57b")) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,  size = 8), 
                axis.title.x = element_blank()) +labs(y="CV")
    ggsave(file.path(output_dir, "boxplot_CV.png"), width=3.5,height=3, dpi=300)
    ggsave(file.path(output_dir, "boxplot_CV.pdf"), width = 3.5, height = 3, dpi = 300)

}

get_CV_4_HKgenes <- function(df, output_dir) {
    dfdata <- df[ 3:length(colnames(df))]
    gene_names <-  df[, "Gene_Name"]
    dfdata <- as.matrix(as.data.frame(dfdata))
    class(dfdata) <- "numeric"
    CV_all <- apply(dfdata, MARGIN = 1, calc_CV)
    mean_all <- apply(dfdata, MARGIN = 1, mean)
    sd_all <- apply(dfdata, MARGIN = 1, sd)
    geomean_all <- apply(dfdata, MARGIN = 1, geoMean)
    alldf <- data.frame(matrix(nrow = length(rownames(dfdata)), ncol = 5))
    colnames(alldf) <- c("mean", "geomean", "sd", "CV", "genes")
    alldf$mean <- as.numeric(mean_all)
    alldf$geomean <- as.numeric(geomean_all)
    alldf$sd <- as.numeric(sd_all)
    alldf$CV <- as.numeric(CV_all)
    alldf$genes <- gene_names
    alldf <- alldf[order(CV_all), ]
    write.csv(alldf, file.path(output_dir, "hkgenesstats.csv"), quote=F, row.names=F)
    return(list("hkstats" = alldf))
}

identifyhkgenes <- function(df, hknum=5, output_dir=getwd()) {
    message("STATUS: identifying ", hknum, " house keeping genes")
    dfdata <- df[df$Class%in% c("Endogenous", "Housekeeping"), 3:length(colnames(df))]
    gene_names <-  df[df$Class %in% c("Endogenous", "Housekeeping"), "Gene_Name"]
    dfdata <-as.matrix(as.data.frame(dfdata))
    class(dfdata) <- "numeric"
    CV_all <- apply(dfdata, MARGIN = 1, calc_CV)
    mean_all <- apply(dfdata, MARGIN = 1, mean)
    sd_all <- apply(dfdata, MARGIN = 1, sd)
    geomean_all <- apply(dfdata, MARGIN = 1, geoMean)
    message("STATUS: Loading data into data.frame")
    alldf <- data.frame(matrix(nrow = length(rownames(dfdata)), ncol = 5))
    colnames(alldf) <- c("mean", "geomean", "sd", "CV", "genes")
    alldf$mean <- as.numeric(mean_all)
    alldf$geomean <- as.numeric(geomean_all)
    alldf$sd <- as.numeric(sd_all)
    alldf$CV <- as.numeric(CV_all)
    alldf$genes <- gene_names
    alldf <- alldf[order(CV_all), ]
    write.csv(alldf, file.path(output_dir, "lowestCVagenesstats.csv"), quote=F, row.names=F)
    message("STATUS: number of genes with lowest CV = ", paste(alldf$genes[1:hknum], collapse=","))
    return(list("hkgenes"= alldf$genes[1:hknum], "allstats"=alldf))
}

violin_plots_hkgenes <- function(countshk, meta.data, group.by="",
    output_dir=getwd(), filename="hkviolinplot", height = 4, width = 5 ) {
    hkmelt <- reshape2::melt(countshk)
    data_variable <- c()
    for (v in hkmelt$variable) {
        data_variable <- c(data_variable, meta.data[meta.data$sample==v, group.by])
    }
    hkmelt$data_variable <- data_variable
    pl <- ggplot(hkmelt, aes(x = Gene_Name, y = value)) +
        geom_violin(alpha = 0.5) +
        geom_jitter(size = 0.5, height = 0, width = 0.1, aes(colour = data_variable)) +
        #   geom_point(hkmelt, aes(color=run), size=0.5, , position = position_jitter(seed = 1, width = 0.2)) +
        theme(legend.position = "none") + labs(y="counts", x="", colour="") +
        theme_Publication() +
        stat_summary(fun = "mean", width = 0.5, geom = "crossbar", color = "red") +
        theme(axis.text.x = element_text(
            angle = 90, vjust = 1, size = 6, hjust = 1
        ))
    ggsave(filename = file.path(output_dir, paste0(filename,".png")), height = height, width = width, dpi = 300)
    return(pl)
}

theme_Publication <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.6, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face = "italic"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

