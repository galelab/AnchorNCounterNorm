#' differential gene expression 
#'
#' Implements a t.test or annova tukey's test to determine DE genes
#' @param exprsdata  normalized log counts
#' @param meta.data meta data for experiment
#' @param compare.column column from meta data you will be comparing (i.e. treatmennt differences the user wants to explore)
#' @param contrastslist list of contrasts of treatments (i.e. c(t"reatment1-control", "treatment2-control"))
#' @param adj.pval.method method used to adjust pvalues for mulitiple hypothesis (default is BH)
#' @param adjbycomparison whether to adjust pvalues by each comparison seperately (default is this) or all together for ttest or wilcox 
#' @param pval.cutoff pvalue cutoff for determining sigificant genes (default is 0.05)
#' @param save.fig whether to save figures in PDF and PNG file format
#' @param check.norm whether to check if data follows  a normal distribution (default is TRUE)
#' @param DE.test which significance test to use (ttest, wilcox or ANOVA/tukeys(default))
#' @param covariate.column column in meta data to use as covariate (ONLY CAN BE USED IN ANOVA TEST)
#' @param output_dir directory to store results 
#' @param volcano.plot whether to make volcano plots showing DE genes (default is TRUE) 
#' @param heatmap.plot whether to make heatmap showing DE genes (default is TRUE)
#' @import stringr
#' @import ggplot2
#' @import ggrepel
#' @import ComplexHeatmap
#' @import circlize
#' 

run_DE_analysis <- function(exprsdata, meta.data, compare.column, contrastslist=NULL, adj.pval.method="BH", adjbycomparison=TRUE,
    pval.cutoff=0.05, save.fig=TRUE, check.norm=TRUE, DE.test="ANOVA", covariate.column=NULL, output_dir=getwd(), volcano.plot=TRUE,
    heatmap.plot=TRUE) {
        genes <- rownames(exprsdata)
        exprsdata <- as.matrix(data.frame(exprsdata, check.names = FALSE))
        class(exprsdata) <- "numeric"
        rownames(exprsdata) <- genes
        if (!(compare.column %in% colnames(meta.data))) {
            stop("WARNING: ", compare.column, " not in header of the meta.data please make sure the name of the column is spelled correctly")
        }
        if (!is.null(covariate.column)) {
            if (!(covariate.column %in% colnames(meta.data))) {
                stop("WARNING: ", covariate.column, " not in header of the meta.data please make sure the name of the column is spelled correctly")
            }
        }
        if (((DE.test=="ttest") || (DE.test=="wilcox")) & is.null(contrastslist)){
            stop("ERROR: if contrastslist aren't specifed then the DE.test parameter\nneeds to be set equal to ANOVA so all comparisons are run ")
        }
        if(isTRUE(check.norm)) {
            S <- shapiro.test(exprsdata)
            if (S$p.value <0.05) {
                message("Data is normally distributed, proceed with T test (1 contrast) or tukeys test (>1 contrast)")
                NORM=TRUE 
            } else {
                message("Data is NOT normally distributed, proceed with wilcox test")
                NORM=FALSE
            }
        }

        tmp <- t(exprsdata)
        totaldata <- merge(tmp, meta.data, by.x="row.names", by.y=1)
        totaldata[, compare.column] <- str_replace(totaldata[, compare.column], "-", "_")
        if ( DE.test=="wilcox"|| DE.test=="ttest" ) {
            message("Running ", DE.test, " test...")
            allresults <- data.frame(matrix(ncol = 7, nrow=length(genes)*length(contrastslist)))
            counter=0
            for (c in contrastslist) {
                message("Running comparison ", c)
                compar <- strsplit(c, "\\s*-\\s*")[[1]]
                for (g in genes) {
                    counter=counter+1
                    gr1 <- totaldata[totaldata[, compare.column] == compar[1], g]
                    gr2 <- totaldata[totaldata[, compare.column] == compar[2], g]
                    if (DE.test=="ttest") {
                        if (isFALSE(NORM) & isTRUE(check.norm)) {
                            message("Data is not normally distributed recommend running wilcox instead")
                        }
                        test <- t.test(gr2, gr1)
                        r <- c(g, c, mean(gr1),mean(gr2), (mean(gr1)-mean(gr2)),
                                test$statistic, test$p.value)
                    } else { 
                        if (isTRUE(NORM)& isTRUE(check.norm)) {
                            message("Data is normally distributed recommend running t test instead")
                        }
                        test <- wilcox.test(gr2, gr1)
                        r <- c(
                            g, c, mean(gr1), mean(gr2), (mean(gr1) - mean(gr2)),
                            test$statistic, test$p.value
                        )
                    }
                    allresults[counter,] <- r
                }
            }
            colnames(allresults)<- c("gene", "comparison", paste0("mean ", compar[1]), 
                paste0("mean ", compar[2]), "LFC", paste0(DE.test, ".statistic"), "pvalue")
            if (adj.pval.method!="none") {
                if (isTRUE(adjbycomparison)) {
                    adjpvals <- c()
                    for (c in contrastslist) {
                        pvals <- allresults[allresults$comparison==c, "pvalue" ]
                        tmpadj <- p.adjust(pvals,method=adj.pval.method, n=length(pvals))
                        adjpvals<-c(adjpvals, tmpadj)
                    }
                    allresults$adj.pval <- adjpvals
                } else {
                    adjpvals <- p.adjust(allresults$pvalue, 
                        method = adj.pval.method, n = length(allresults$pvalue))
                    allresults$adj.pval <- adjpvals
                }
                allresultssig <- allresults[allresults$adj.pval<=pval.cutoff, ]
            } else {
                allresultssig <- allresults[allresults$pvalue<=pval.cutoff, ]
            }
            write.csv(allresults, file.path(output_dir, paste0(DE.test, "_allResults.csv")))
            write.csv(allresultssig, file.path(output_dir, paste0(DE.test, "_SignificantResults.csv")))
            allresults4fig <- allresults[, c("gene","LFC","comparison","adj.pval" )]
            allresults4figsig <- allresultssig[, c("gene","LFC","comparison","adj.pval" )]
            if(nrow(allresults4figsig) > 0 ) {
                if (isTRUE(volcano.plot)) {
                    volcanoplotfigures <- volcano_plot(allresults4fig, pval.cutoff = pval.cutoff, save.fig = save.fig, output_dir = output_dir)
                }
                if ((isTRUE(heatmap.plot) & length(unique(allresults4figsig$comparison))>1)) {
                    heatmapfigures <- heatmap4DE(allresults4fig, allresults4figsig, contrastslist,
                        save.fig=save.fig, output_dir = output_dir)
                } else if ((isTRUE(heatmap.plot) & length(unique(allresults4figsig$comparison))==1)) {
                    message('WARNING: heatmap will not be made with only one comparison')
                }
            } else {
                message("WARNING: no significant genes ")
                if (isTRUE(volcano.plot)) {
                    volcanoplotfigures <- volcano_plot(allresults4fig,
                        pval.cutoff = pval.cutoff,
                        save.fig = save.fig, output_dir = output_dir
                    )
                }
            }
            if(nrow(allresults4figsig) > 0 ) {
                if (isTRUE(volcano.plot) & (isTRUE(heatmap.plot))) {
                    return(list("allresults" <- allresults4fig, "sigresults" = allresults4figsig, "volcanoplots" = volcanoplotfigures, "heatmapplots" = heatmapfigures))
                } else if ((isTRUE(volcano.plot) & isFALSE(heatmap.plot))) {
                    return(list("allresults" <- allresults4fig, "sigresults" = allresults4figsig, "volcanoplots" = volcanoplotfigures))
                } else if ((isFALSE(volcano.plot) & isTRUE(heatmap.plot))) {
                    return(list("allresults" <- v, "sigresults"=allresults4figsig, "heatmapplots"=heatmapfigures))
                } else {
                    return(list("allresults" <- allresults4fig, "sigresults" = allresults4figsig))
                }
            } else {
                return(list("allresults" <- allresults4fig, "sigresults" = allresults4figsig))
            }
        } else if (DE.test=="ANOVA") {
            comparisonlist <- list()
            totaldatamelt <- reshape2::melt(totaldata)
            totaldatamelt$value <- as.numeric( totaldatamelt$value)
            for (g in genes ){
                totaldatamelttmp <- totaldatamelt[totaldatamelt$variable==g, ]
                if (is.null(covariate.column)) {
                    model=lm( totaldatamelttmp[, "value"] ~ totaldatamelttmp[, compare.column] )
                } else {
                    model=lm( totaldatamelttmp[, "value"] ~ totaldatamelttmp[, compare.column] + totaldatamelttmp[, covariate.column] )
                }
                ANOVA=aov(model)
                TUKEY <- TukeyHSD(x=ANOVA, 'totaldatamelttmp[, compare.column]', conf.level=0.95)
                T <- as.data.frame(TUKEY[[1]])
                T$comparison <- rownames(T)
                T$gene <- rep(g, nrow(T))
                comparisonlist[[g]] <- as.data.frame(T)
            }
            alltukeydf <- do.call(rbind, comparisonlist)
            rownames(alltukeydf) <- 1:nrow(alltukeydf)
            if (!is.null(contrastslist)) {
                contrastslist <- str_replace_all(contrastslist, "\\s*-\\s*", "-")
                alltukeydf<- alltukeydf[alltukeydf$comparison %in% contrastslist, ]
            }
            colnames(alltukeydf)[4] <- "adj.pval"
            colnames(alltukeydf)[1] <- "LFC"
            alltukeydfsig <- alltukeydf[alltukeydf[, "adj.pval"]<=pval.cutoff, ]
            write.csv(alltukeydf, file.path(output_dir, "TukeysallResults.csv"))
            write.csv(alltukeydfsig, file.path(output_dir, "TukeysSignificantResults.csv"))
            alltukeydf4fig <- alltukeydf[, c("gene","LFC","comparison","adj.pval" )]
            if(nrow(alltukeydfsig) > 0 ) {
                if (isTRUE(volcano.plot)) {
                    volcanoplotfigures <- volcano_plot(alltukeydf4fig, pval.cutoff = pval.cutoff, 
                        save.fig = save.fig, output_dir = output_dir)
                }
                if ((isTRUE(heatmap.plot) & length(unique(alltukeydfsig$comparison))>1)) {
                    heatmapfigures <- heatmap4DE(alltukeydf4fig, alltukeydfsig, contrastslist,
                        save.fig = save.fig, output_dir = output_dir)
                } else if ((isTRUE(heatmap.plot) & length(unique(alltukeydfsig$comparison))==1)) {
                    message('WARNING: heatmap will not be made with only one comparison')
                }
            } else {
                message("WARNING: no significant genes ")
                if (isTRUE(volcano.plot)) {
                    volcanoplotfigures <- volcano_plot(alltukeydf4fig,
                        pval.cutoff = pval.cutoff,
                        save.fig = save.fig, output_dir = output_dir
                    )
                }
            }
            if(nrow(alltukeydfsig) > 0 ) {
                if (isTRUE(volcano.plot) & (isTRUE(heatmap.plot)) & length(unique(alltukeydfsig$comparison))>1) {
                    return(list("allresults" = alltukeydf, "sigresults" = alltukeydfsig, "volcanoplots" = volcanoplotfigures, "heatmapplots" = heatmapfigures))
                } else if ((isTRUE(volcano.plot) & isFALSE(heatmap.plot))) {
                    return(list("allresults" = alltukeydf, "sigresults" = alltukeydfsig, "volcanoplots" = volcanoplotfigures))
                } else if ((isFALSE(volcano.plot) & comparisonlistcomparisonlist(heatmap.plot) & length(unique(alltukeydfsig$comparison))>1)) {
                    return(list("allresults" = alltukeydf, "sigresults" = alltukeydfsig, "heatmapplots" = heatmapfigures))
                } else {
                    return(list("allresults" = alltukeydf, "sigresults" = alltukeydfsig))
                }
            } else {
                return(list("allresults" = alltukeydf, "sigresults" = alltukeydfsig))
            }
        } else {  
            stop("ERROR: DE.test needs to either be ttest, wilcox or ANOVA")
        }
    }

volcano_plot <- function(df, pval.cutoff, save.fig=TRUE, output_dir=getwd()) {
    df <- as.data.frame(df)
    individual.compare.plots <- list()
    for (comp in unique(df$comparison)) {
        tmp <- df[df$comparison==comp, ]
        tmp$Sign <- "NS"
        tmp$Sign[tmp$adj.pval < pval.cutoff] <- paste0("Padj < ", pval.cutoff)
        tmp$Sign[tmp$adj.pval < 0.01] <- "Padj < 0.01"
        significance <- paste0("Padj < ", pval.cutoff)
        tmp$Sign <- factor(tmp$Sign,
            levels = c(
                "NS", paste0("Padj < ", pval.cutoff), "Padj < 0.01"
            )
        )
        tmp$LFC <- as.numeric( tmp$LFC)
        tmp$adj.pval <- as.numeric( tmp$adj.pval)
        pl <- ggplot(tmp, aes(x=LFC, y= -log10(adj.pval),  color = Sign, label = gene)) + 
            geom_vline(xintercept = 0, lty = "dashed",  linewidth=0.2) +
            geom_hline(yintercept = -log10(pval.cutoff), lty = "dashed",  linewidth=0.2) +
            geom_point(size = 3, alpha=0.5) +
            labs(
                x = "log2(FC)",
                y = "Significance, -log10(P)",
                color = "Sig"
            ) +
            scale_color_manual(
                values = c(
                   "gray","dodgerblue","orange"
                ),
                guide = guide_legend(override.aes = list(size = 2))
            ) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
            geom_text_repel(
                data =subset(tmp, adj.pval  < as.numeric(pval.cutoff)), aes(label=gene),
                size = 2,  point.padding = 0.15, color = "black",
                min.segment.length = .1, box.padding = .2, lwd = 1,
                max.overlaps = 10
            ) +
            theme_minimal(base_size = 12) +
            theme(legend.title =element_blank(), legend.position = "bottom",legend.key.size = unit(0.4, "cm"), 
                axis.line = element_line(colour = "black", linewidth=0.2),
                panel.grid.major = element_line(colour = "#f0f0f0"),
                panel.grid.minor = element_blank()) +
                theme(strip.background = element_rect(fill = "white")) 
        if(isTRUE(save.fig)) {
            ggsave(file.path(output_dir, paste0("volcanoplot_individual_", comp,".png")), width=3, height=3.5, units="in", dpi = 300, bg = "white")
            ggsave(file.path(output_dir, paste0("volcanoplot_individual_", comp,".pdf")), width=3, height=3.5, units="in",  dpi = 300, bg = "white")
        }
    }
    if(length(unique(df$comparison)) > 1) {
        df$Sign <- "NS"
        df$Sign[df$adj.pval < pval.cutoff] <- paste0("Padj < ", pval.cutoff)
        df$Sign[df$adj.pval < 0.01] <- "Padj < 0.01"
        significance <- paste0("Padj < ", pval.cutoff)

        df$Sign <- factor(df$Sign,
            levels = c(
                "NS", paste0("Padj < ", pval.cutoff), "Padj < 0.01"
            ))
        df$LFC <- as.numeric( df$LFC)
        df$adj.pval <- as.numeric(df$adj.pval)

        pl <- ggplot(df, aes(x=LFC, y=-log10(adj.pval),  color = Sign, label = gene)) +
            geom_vline(xintercept = 0, lty = "dashed",  linewidth=0.2) +
            geom_hline(yintercept = -log10(pval.cutoff), lty = "dashed",  linewidth=0.2) +
            geom_point(size = 3, alpha=0.5) +
            facet_wrap(~comparison, ncol=4, labeller = labeller(.default = label_both)) +
            labs(
                x = "log2(FC)",
                y = "Significance, -log10(P)",
                color = "Sig"
            ) +
            scale_color_manual(
                values = c(
                   "gray","dodgerblue","orange"

                ),
                guide = guide_legend(override.aes = list(size = 2))
            ) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
            geom_text_repel(
                data =subset(df, adj.pval  < as.numeric(pval.cutoff)), aes(label=gene),
                size = 2,  point.padding = 0.15, color = "black",
                min.segment.length = .1, box.padding = .2, lwd = 1,
                max.overlaps = 10
            ) +
            theme_minimal(base_size = 12) +
            theme(legend.title =element_blank(), legend.position = "bottom",legend.key.size = unit(0.4, "cm"), 
                axis.line = element_line(colour = "black", linewidth=0.2),
                panel.grid.major = element_line(colour = "#f0f0f0"),
                panel.grid.minor = element_blank()) +
                theme(strip.background = element_rect(fill = "white")) 
        if(isTRUE(save.fig)) {
            ggsave(file.path(output_dir, paste0("volcanoplot_all_", comp,".png")), dpi = 300, bg = "white")
            ggsave(file.path(output_dir, paste0("volcanoplot_all_", comp,".pdf")), dpi = 300, bg = "white")
        }
    }
     return(list("individual.plots"=individual.compare.plots, "combined.comparison.plots"=pl))
}

heatmap4DE <- function(df, dfsig, contrastsorder, save.fig=TRUE, output_dir=getwd()) {
    dfcom <- reshape2::dcast(data = df, formula = gene~comparison, value.var = "LFC")
    dfcomsigonly <- reshape2::dcast(data = dfsig, formula = gene~comparison, value.var = "LFC")
    # Sig data frame 
    dfcomsigonly4sig <- dfcomsigonly
    rownames(dfcomsigonly4sig) <- dfcomsigonly4sig$gene
    dfcomsigonly4sig$gene <- NULL
    write.csv(dfcomsigonly4sig, file.path(output_dir, "DEgenesSigLFConly.csv"))
    dfcomsigonly4sig <- as.matrix(dfcomsigonly4sig)
    dfcomsigonly4sig <- dfcomsigonly4sig[, intersect(contrastsorder, colnames(dfcomsigonly4sig))]
    class(dfcomsigonly4sig) <- "numeric"
    # Full data frame 
    dfcomsig4hm <- dfcom[dfcom$gene %in% dfcomsigonly$gene, ]
    rownames(dfcomsig4hm) <- dfcomsig4hm$gene
    dfcomsig4hm$gene <- NULL
    write.csv(dfcomsig4hm, file.path(output_dir, "DEgenesSigLFC.csv"))
    dfcomsig4hm <- as.matrix(dfcomsig4hm)
    dfcomsig4hm <- dfcomsig4hm[, contrastsorder]
    class(dfcomsig4hm) <- "numeric"

    col_fun <- colorRamp2(c(-3, -2, -1, 0, 1, 2, 3), c("darkblue", "mediumblue", "dodgerblue", "white", "orange", "red", "darkred"))
    hmC <- Heatmap(dfcomsig4hm,
        name = "LFC",
        col = col_fun,
        show_row_names = TRUE,
        cluster_row_slices = T,
        cluster_columns = FALSE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        cluster_column_slices = FALSE,
        column_names_rot = 90, 
        # row_labels = rownames(dfcomsig4hm),
        border = F, width = unit(4.5, "cm"),
        height = unit(4.5, "cm"),
        column_names_gp = grid::gpar(fontsize = 7),
        column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
        row_names_gp = grid::gpar(fontsize = 6),
        heatmap_legend_param = list(
            legend_height = unit(3, "cm"),
            labels_gp = gpar(fontsize = 14),
            title_gp = gpar(fontsize = 14, fontface = "bold")
        ), use_raster = TRUE, raster_quality = 5)

    dfcomsigonly4sig[is.na(dfcomsigonly4sig)] <- 0
    hmCsig <- Heatmap(dfcomsigonly4sig,
        name = "LFC",
        col = col_fun,
        show_row_names = TRUE,
        cluster_row_slices = T,
        cluster_columns = FALSE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        cluster_column_slices = FALSE,
        column_names_rot = 90, 
        # row_labels = rownames(dfcomsig4hm),
        border = F, width = unit(4.5, "cm"),
        height = unit(4.5, "cm"),
        column_names_gp = grid::gpar(fontsize = 7),
        column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
        row_names_gp = grid::gpar(fontsize = 6),
        heatmap_legend_param = list(
            legend_height = unit(3, "cm"),
            labels_gp = gpar(fontsize = 14),
            title_gp = gpar(fontsize = 14, fontface = "bold")
        ), use_raster = TRUE, raster_quality = 5)

        if (isTRUE(save.fig)) {
            pdf(file.path(output_dir, "Heatmap.pdf")) # , width = 6, height = 6, units = "in", res = 700)
            heatmapCOV <- draw(hmC,
                heatmap_legend_side = "left"
            )
            dev.off()
            png(file.path(output_dir, "Heatmap.png"), width = 6, height = 6, units = "in", res = 700)
            heatmapCOV <- draw(hmC,
                heatmap_legend_side = "left"
            )
            dev.off()
            pdf(file.path(output_dir, "HeatmapSigLFCs.pdf")) # , width = 6, height = 6, units = "in", res = 700)
            heatmapCOVsig <- draw(hmCsig,
                heatmap_legend_side = "left"
            )
            dev.off()
            png(file.path(output_dir, "HeatmapSigLFCs.png"), width = 6, height = 6, units = "in", res = 700)
            heatmapCOVsig <- draw(hmCsig,
                heatmap_legend_side = "left"
            )
            dev.off()
        } else {
            heatmapCOV <- draw(hmC,
                heatmap_legend_side = "left"
            )
            heatmapCOVsig <- draw(hmCsig,
                heatmap_legend_side = "left"
            )
        }
    return(list("hmsiggenes"=heatmapCOV, "hmsiggenes.siglfcsonly"=heatmapCOVsig))
}