#' dimentionality reduction 
#'
#' This implements pca or umap dimensionality reduction on data from normalization step 
#' @param counts  matrix counts
#' @keywords house keeping genes
#' @export
#' @import factoextra
#' @import Rtsne
#' @import ggplot2
#' @import Polychrome
#' @examples
#' data <- hk_gene_stats(counts)
#' 
#' 

dim_reduction <- function(exprs, meta.data,  target_columns=c(1,2),
    reduction="PCA", save.fig=TRUE, output_dir=getwd(), pointsize=3) {
    c1ass1 = meta.data[, target_columns[1]]
    class2 = meta.data[, target_columns[2]]
    P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))

    if (reduction=="PCA") {
        pca <- prcomp(t(exprs))
        E <- get_eig(pca)
        cx <- sweep(t(exprs), 2, colMeans(t(exprs)), "-")
        PCA <- pca$x
        minx <- min(PCA[, 1])
        maxx <- max(PCA[, 1])
        miny <- min(PCA[, 2])
        maxy <- max(PCA[, 2])

        pcaplot <- ggplot(aes(x=PCA[, 1], y=PCA[, 2], color = factor(class2), shape = factor(class1), size = I(pointsize))) +
            theme_minimal() +
            theme(legend.title = element_blank()) +
            xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
            ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
            theme(legend.position = "right") +
            scale_color_manual(values = as.character(P36)) +
            scale_fill_manual(values = as.character(P36)) +
            scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        if (isTRUE(save.fig)) {
            ggsave(ggsave(file.path(output_dir, paste0(reduction,".png")), width = 4.5, height = 3, bg = "white", dpi = 300))
            ggsave(ggsave(file.path(output_dir, paste0(reduction, ".pdf")), width = 4.5, height = 3, bg = "white", dpi = 300))
        }
        #scree plot 
        scree.plot <- fviz_eig(PCA, addlabels = TRUE, hjust = -0.3)
        if (isTRUE(save.fig)) {
            png(file.path(output_dir, paste0(reduction, "scree.png")), width = 7, height = 6, units = "in", res = 300)
            print(scree.plot)
            dev.off()
            pdf(file.path(output_dir, paste0(reduction, "scree.pdf")), width = 7, height = 6, units = "in", res = 300)
            print(scree.plot)
            dev.off()
        }

        loadingscores <- as.data.frame(pca$rotation)
        is_pc1_0 <- loadingscores$PC1 > 0
        is_pc2_0 <- loadingscores$PC2 > 0

        loadingscores <- loadingscores[is_pc1_0, ]
        loadingscores <- loadingscores[with(loadingscores, order(-PC1)), ]

        write.table(loadingscores["PC1"], file = file.path(output_dir, "loadingscores_pc1.txt"))

        loadingscores <- as.data.frame(pca$rotation)
        loadingscores <- loadingscores[is_pc2_0, ]
        loadingscores <- loadingscores[with(loadingscores, order(-PC2)), ]
        write.table(loadingscores["PC2"], file = file.path(output_dir, "loadingscores_pc2.txt"))
    } else if (reduction=="UMAP") {
        UMAP <- umap(t(exprs))
        U <- umap$layout
        minx <- min(U[, 1])
        maxx <- max(U[, 1])
        miny <- min(U[, 2])
        maxy <- max(U[, 2])

        umapplot <- ggplot(aes(x = U[, 1], y = U[, 2], color = factor(class2), shape = factor(class1), size = I(pointsize))) +
            theme_minimal() +
            theme(legend.title = element_blank()) +
            theme(legend.position = "right") +
            scale_color_manual(values = as.character(P36)) +
            scale_fill_manual(values = as.character(P36)) +
            scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        if (isTRUE(save.fig)) {
            ggsave(ggsave(file.path(output_dir, paste0(reduction, ".png")), width = 4.5, height = 3, bg = "white", dpi = 300))
            ggsave(ggsave(file.path(output_dir, paste0(reduction, ".pdf")), width = 4.5, height = 3, bg = "white", dpi = 300))
        }
    } else {
        stop("ERROR: reduction parameter has to set to PCA or UMAP")
    }
    if (reduction=="PCA") {
        return(list("pcaobject" = pca, "pcafigure"=pcaplot, "screeplot"=scree.plot))
    } else if (reduction=="UMAP") {
        return(list("umapobject"=UMAP, "umapfigure"=umapplot))
    }
}