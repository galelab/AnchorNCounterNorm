#' dimentionality reduction 
#'
#' This implements pca or umap dimensionality reduction on data from normalization step 
#' @param exprsdata  normalized log counts
#' @param meta.data meta data for experiment
#' @param target_columns meta data to display by shape and color (2 columns can be shown)
#' @param reduction what feature reduction algorithm to use PCA (default) or UMAP
#' @param save.fig whether to save figures in a file (defualt is TRUE)
#' @param output_dir where to save results (default is current working directory)
#' @param pointsize size of dots that represent samples on pca or umap (default is 3)
#' @keywords house keeping genes
#' @export
#' @import factoextra
#' @import umap
#' @import ggplot2
#' @import Polychrome
#' @import grDevices
#' @examples
#' data <- hk_gene_stats(counts)
#' 
#' 

dim_reduction <- function(exprsdata, meta.data,  target_columns=c(2,3),
    reduction="PCA", save.fig=TRUE, output_dir=getwd(), file.name="PCA", pointsize=3, 
    ordertargetcolumn1=NULL, ordertargetcolumn2=NULL, colorgradient=FALSE) {

    class1 = meta.data[, target_columns[1]]
    class2 = meta.data[, target_columns[2]]
    exprsdata <- as.matrix(as.data.frame(exprsdata))
    class(exprsdata) <-"numeric"
    if (reduction=="PCA") {
        pca <- prcomp(t(exprsdata))
        E <- get_eig(pca)
        PCA <- pca$x
        df <- as.data.frame(cbind(PCA[, 1], PCA[, 2], as.character(class1), as.character(class2)))
        colnames(df) <- c("PC1", "PC2", "class1", "class2")
        df$PC1 <- as.numeric(df$PC1)
        df$PC2 <- as.numeric(df$PC2)
        if (!is.null(ordertargetcolumn1)) {
            df$class1 <- factor(df$class1, levels = ordertargetcolumn1)
        } else {
            df$class1 <- factor(df$class1)
        }
        if (!is.null(ordertargetcolumn2)) {
            df$class2 <- factor(df$class2, levels = ordertargetcolumn2)
        } else {
            df$class2 <- factor(df$class2)
        }
        if (isFALSE(colorgradient)) {
            P36 <- createPalette(length(levels(df$class2)), c("#ff0000", "#00ff00", "#0000ff"))
        } else if( isTRUE(colorgradient)) {
            colfunc <- colorRampPalette(c("#bcb7b7","#797373","#00ff00", "#015601"))
            P36 <- colfunc(length(levels(df$class2)))
        }
        pcaplot <- ggplot(df, aes(x=PC1, y=PC2, color = class2, shape = class1)) +
            geom_point(size = I(pointsize)) +
            theme_minimal() +
            theme(legend.title = element_blank()) +
            xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
            ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
            theme(legend.position = "right") + theme(legend.key.size = unit(0.2, "cm"),
                panel.grid.major = element_line(colour = "#f0f0f0"),
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black", linewidth=0.2)) +
            scale_color_manual(values = as.character(P36)) +
            scale_fill_manual(values = as.character(P36)) +
            scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        if (isTRUE(save.fig)) {
            ggsave(file.path(output_dir, paste0(file.name, "_", reduction,".png")), width = 4.5, height = 3, bg = "white", dpi = 300)
            ggsave(file.path(output_dir, paste0(file.name, "_", reduction, ".pdf")), width = 4.5, height = 3, bg = "white", dpi = 300)
        }
        #scree plot 
        scree.plot <- fviz_eig(pca, addlabels = TRUE, hjust = -0.3)
        if (isTRUE(save.fig)) {
            png(file.path(output_dir, paste0(file.name, "_", reduction, "scree.png")), width = 7, height = 6, units = "in", res = 300)
            print(scree.plot)
            dev.off()
            pdf(file.path(output_dir, paste0(file.name, "_", reduction, "scree.pdf")), width = 7, height = 6)
            print(scree.plot)
            dev.off()
        }

        loadingscores <- as.data.frame(pca$rotation)
        is_pc1_0 <- loadingscores$PC1 > 0
        is_pc2_0 <- loadingscores$PC2 > 0

        loadingscores <- loadingscores[is_pc1_0, ]
        loadingscores <- loadingscores[with(loadingscores, order(-PC1)), ]

        write.table(loadingscores["PC1"], file = file.path(output_dir, paste0(file.name, "_","loadingscores_pc1.txt")))

        loadingscores <- as.data.frame(pca$rotation)
        loadingscores <- loadingscores[is_pc2_0, ]
        loadingscores <- loadingscores[with(loadingscores, order(-PC2)), ]
        write.table(loadingscores["PC2"], file = file.path(output_dir, paste0(file.name, "_","loadingscores_pc2.txt")))
    } else if (reduction=="UMAP") {
        UMAP <- umap(t(exprsdata), n_neighbors=5)
        U <- UMAP$layout
        df <- as.data.frame(cbind(U[, 1], U[, 2], class1, class2))
        colnames(df) <- c("UMAP1", "UMAP2", "class1", "class2")
        df$UMAP1 <- as.numeric(df$UMAP1)
        df$UMAP2 <- as.numeric(df$UMAP2)
        umapplot <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = factor(class2), 
            shape = factor(class1))) +
            geom_point(size=pointsize) +
            theme_minimal() +
            theme(legend.title = element_blank()) +
            theme(legend.position = "right") + theme(legend.key.size = unit(0.2, "cm"),
                panel.grid.major = element_line(colour = "#f0f0f0"),
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black", linewidth=0.2)) + 
            scale_color_manual(values = as.character(P36)) +
            scale_fill_manual(values = as.character(P36)) +
            scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        if (isTRUE(save.fig)) {
            ggsave(file.path(output_dir, paste0(file.name, "_", reduction, ".png")), width = 4.5, height = 3, bg = "white", dpi = 300)
            ggsave(file.path(output_dir, paste0(file.name, "_", reduction, ".pdf")), width = 4.5, height = 3, bg = "white", dpi = 300)
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