#' Load nCounter files 
#'
#' This loads count and target file info (internal function only)
#' @param pathtoRCC  matrix raw or normalized (default is current working directory)
#' @param target_file file with experimental information (meta data csv file)
#' @keywords Ncounter RCC files
#' @export
#' @import NanoStringQCPro
#' @import ggplot2
#' @import reshape2
#' @examples
#' data <- load_nCounter_files(pathtoRCC="./", meta.data="targetfile.csv")
#' 
load_nCounter_files <- function(pathtoRCC=".", meta.data="", save.fig=TRUE) {
    # read in meta data 
    target = read.csv(meta.data, header = TRUE)
    if (!("sample" %in% colnames(target))) {
        message("WARNING: if header doesn't contain column named, sample, please rename in meta data file or specify column with sample ids in future steps with sample_column option. sample column should have rcc file names")
    }

    # read in rcc files 
    files.RCC = list.files(path=pathtoRCC, pattern="RCC")

    # check to ensure RCC files are present 
    if (length(files.RCC) > 0 ) { 
        message("STATUS: Number of RCC files is ", length(files.RCC))
    } else {
        stop("No rcc files found in directory ", pathtoRCC, " fix before continuing")
    }

    # generate data frames to store data
    raw_counts = as.data.frame(matrix(nrow = 0, ncol = length(files.RCC) + 2))
    colnames(raw_counts)[1:2] = c("Gene_Name", "Class")
    qc_nCounter_data = as.data.frame(matrix(nrow = length(files.RCC), ncol = 11))
    colnames(qc_nCounter_data) = c(
        "BCAC_ID", "SampleID", "Owner", "Comments", 
        "Date", "GeneRLF", "SystemAPF", "imagingQC",
        "bindingDensityQC", "limitOfDetectionQC", "positiveLinearityQC"
    )

    # Load in RCC files and put counts in raw_counts (dataframe) and qc info in qc_nCounter_data (dataframe)
    for (i in 1:length(files.RCC)) {
        message("STATUS: loading file ", files.RCC[i])
        rcc = readRcc(file.path(pathtoRCC, files.RCC[i]))

        if (i==1) {
            #When loading the first RCC file generate rows for each gene in panel
            raw_counts[nrow(raw_counts) + length(rcc$Code_Summary$Count), ] <- NA
            
            #fill Gene and class columns in data frame with gene and class info from the first RCC file
            raw_counts$Gene_Name <- as.character(rcc$Code_Summary$Name)
            raw_counts$Class <- as.character(rcc$Code_Summary$CodeClass)
        } else if (length(rcc$Code_Summary$Count) != dim(raw_counts)[1]){
            #This checks to make sure subsequent RCC files have the same number of genes as first one
            stop("Number of genes in RCC file ", length(rcc$Code_Summary$Count), 
                "does not match previous file ", dim(raw_counts)[1], 
                ". RCC files need to have the same number of genes if being merged together")

        } else if (isFALSE(all.equal(raw_counts$Gene_Name, as.character(rcc$Code_Summary$Name)))) {
            #This checks to make sure genes in subsequent RCC file are in the same order if not tries to change 
            #order of genes to the same as the first RCC file. If code can't it will error out
            message("STATUS: Genes are not in the same order as previous RCC file attempting to fix...")
            rcc$Code_Summary = rcc$Code_Summary[match(raw_counts$Gene_Name, rcc$Code_Summary$Name), ]

            if (isTRUE(all.equal(raw_counts$Gene_Name, as.character(rcc$Code_Summary$Name)))) {
                message("STATUS: fixed order continuing... ")
            } else {
                stop("Could not fix gene order issue please fix manually in rcc files")
            }
        }

        # Load data from rcc file into data frames if no errors have arisen
        raw_counts[, i + 2] <- as.numeric(rcc$Code_Summary$Count)
        colnames(raw_counts)[i + 2] <- files.RCC[i]
        qc_nCounter_data[i, 2:7] <- as.vector(rcc$Sample_Attributes)
        qc_nCounter_data$imagingQC[i] <- imagingQC(rcc)
        qc_nCounter_data$bindingDensityQC[i] <- bindingDensityQC(rcc, .05, 2.25)
        qc_nCounter_data$limitOfDetectionQC[i] <- limitOfDetectionQC(rcc)
        qc_nCounter_data$positiveLinearityQC[i] <- positiveLinQC(rcc)
        rownames(qc_nCounter_data)[i] <- files.RCC[i]
    }

    # Generate figures to examine counts
    melted_raw <- melt(raw_counts)
    ggplot(melted_raw, aes(x = variable, y = value, fill = Class)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        scale_fill_manual(values = c("#ff615d", "#004c7a", "#ffd400", "#61c57b")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
        ggtitle("Total raw counts") + labs(x="", y="raw counts")
    if(isTRUE(save.fig)) {
        message("STATUS: saving figure as png file here: ", getwd())
        ggsave(file.path(getwd(), "barplot_totalrawcounts.png"), bg="white", dpi=300)
    }

    ggplot(melted_raw, aes(x = variable, y = log2(value), fill = Class)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        scale_fill_manual(values = c("#ff615d", "#004c7a", "#ffd400", "#61c57b")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,  size = 8)) +
        ggtitle("Log2() raw counts") + labs(x = "", y = "log2 raw counts") 

    if (isTRUE(save.fig)) {
        message("STATUS: saving figure as png file here: ", getwd())
        ggsave(file.path(getwd(), "barplot_log2rawcounts.png"), bg = "white", dpi = 300)
    }

    # raw counts all sample counts - faceted barplot by gene class
    ggplot(melted_raw, aes(x = variable, y = value, fill = Class)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        scale_fill_manual(values = c("#ff615d", "#004c7a", "#ffd400", "#61c57b")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,  size = 8)) +
        ggtitle("Total raw counts") + labs(x = "", y = "raw counts") +
        facet_wrap(~Class, nrow = 1)
    if (isTRUE(save.fig)) {
        message("STATUS: saving figure as png file here: ", getwd())
        ggsave(file.path(getwd(), "barplot_Class_totalrawcounts.png"), bg = "white", dpi = 300)
    }
    # log2() all sample counts - faceted barplot by gene class
    ggplot(melted_raw, aes(x = variable, y = log2(value), fill = Class)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        scale_fill_manual(values = c("#ff615d", "#004c7a", "#ffd400", "#61c57b")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,  size = 8)) +
        ggtitle("e) Log2() raw counts") + labs(x = "", y = "log2 raw counts") +
        facet_wrap(~Class, nrow = 1)
    if (isTRUE(save.fig)) {
        message("STATUS: saving figure as png file here: ", getwd())
        ggsave(file.path(getwd(), "barplot_Class_log2rawcounts.png"), bg = "white", dpi = 300)
    }
    write.csv(raw_counts, file.path(getwd(), "counts.csv"))
    write.csv(raw_counts, file.path(getwd(), "qc.csv"))

    return(list(counts = raw_counts, meta.data = target, qc = qc_nCounter_data))

}
