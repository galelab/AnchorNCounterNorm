# revolutionaryGameChanger
Nanostring anchoring normalization package


## Install 
To download in R please use install_github('galelab/revolutionaryGameChanger').  First load library('devtools')

Then to load library library(NCounterNorm)

## How to execute pipeline

data = load_nCounter_files(pathtoRCC="folder/with/RCCfiles/", meta.data="path/to/metadata/csvfile")

HKstats <- hk_gene_stats(data, manually_selected_hk_genes = c("B2M", "GAPDH", "PGK1", "RPLP0"), group.by="day")

normdata <- ratio_normalization(data, hkgenes=c("GAPDH", "PGK1", "RPLP0"))

dimred <- dim_reduction(normdata$log_counts_ratio, data$meta.data, target_columns = c(2,4))

### DE examples 

DE1 <- run_DE_analysis(normdata$log_counts_ratio, data$meta.data, compare.column = "day", 
    contrastslist=c("day_3-day_0"), DE.test="ttest")

DE2 <- run_DE_analysis(normdata$log_counts_ratio, data$meta.data, compare.column = "animalID", pval.cutoff = 0.05, contrastslist=c("animal2-animal1", "animal3-animal1"), DE.test="ANOVA")