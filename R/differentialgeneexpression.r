#' differential gene expression 
#'
#' Implements a t.test or annova tukey's test to determine DE genes
#' @param exprsdata  normalized log counts
#' @param meta.data meta data for experiment
#' 
#' 

run_DE_analysis <- function(exprsdata, meta.data, meta.data.column, contrasts=c(), adj.pval.method="BH", adjbycomparison=TRUE,
    pval.cutoff=0.05, save.fig=TRUE, check.norm=TRUE, DE.test=NULL, output_dir=getwd()) {
        genes <- rownames(exprsdata)
        exprsdata <- as.matrix(data.frame(exprsdata, check.names = FALSE))
        class(exprsdata) <- "numeric"
        rownames(exprsdata) <- genes
        if(isTRUE(check.norm) & is.null(DE.test)) {
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

        if ( DE.test=="wilcox"|| DE.test=="ttest" ) {
            message("Running ", DE.test, " test...")
            allresults <- data.frame(matrix(ncol = 7, nrow=length(genes)*length(contrasts)))
            counter=0
            for (c in contrasts) {
                message("Running comparison ", c)
                compar <- strsplit(c, "\\s*-\\s*")[[1]]
                for (g in genes) {
                    counter=counter+1
                    gr1 <- totaldata[totaldata[,meta.data.column]==compar[1], g]
                    gr2 <- totaldata[totaldata[, meta.data.column] == compar[2], g]
                    if (DE.test=="ttest") {
                        if (isFALSE(NORM)) {
                            message("Data is not normally distributed recommend running wilcox instead")
                        }
                        test <- t.test(gr2, gr1)
                        r <- c(g, c, mean(gr1),mean(gr2), (mean(gr1)-mean(gr2)),
                                test$statistic, test$p.value)
                    } else { 
                        if (isTRUE(NORM)) {
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
                    for (c in contrast) {
                        pvals <- allresults[allresults$comparisons==c, "pvalue" ]
                        tmpadj <- p.adjust(pvals,method=adj.pval.method, n=length(pvals))
                        adjpvals<-c(adjpvals, tmpadj)
                    }
                    allresults$adj.pval <- adjpvals
                } else {
                    adjpvals <- p.adjust(allresults$pvalue, 
                        method = adj.pval.method, n = length(allresults$pvalue))
                    allresults$adj.pval <- adjpvals
                }
            }
        }
    }