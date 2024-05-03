# AnchorNCounterNorm
NCounter anchoring normalization package

## Install 
To download in R please use install_github('galelab/AnchorNCounterNorm').  First load library('devtools')

```R
library('devtools')
install_github('galelab/AnchorNCounterNorm')
library(AnchorNCounterNorm)
```
## How to execute pipeline

```R
data <- load_nCounter_files(pathtoRCC="folder/with/RCCfiles/", meta.data="path/to/metadata/csvfile")

HKstats <- hk_gene_stats(data, manually_selected_hk_genes = c("B2M", "GAPDH", "PGK1", "RPLP0"), group.by="day")

normdata <- ratio_normalization(data, hkgenes=c("GAPDH", "PGK1", "RPLP0"))

dimred <- dim_reduction(normdata$log_counts_ratio, data$meta.data, target_columns = c(2,4))
```
### DE examples 

```R
DE1 <- run_DE_analysis(normdata$log_counts_ratio, data$meta.data, compare.column = "day", 
    contrastslist=c("day_3-day_0"), DE.test="ttest")

DE2 <- run_DE_analysis(normdata$log_counts_ratio, data$meta.data, compare.column = "animalID", pval.cutoff = 0.05, contrastslist=c("animal2-animal1", "animal3-animal1"), DE.test="ANOVA")
```
## Citation 
Method was described in the following publication:

Davis MA, Voss K, Turnbull JB, Gustin AT, Knoll M, Muruato A, Hsiang TY, Dinnon Iii KH, Leist SR, Nickel K, Baric RS, Ladiges W, Akilesh S, Smith KD, Gale M Jr. A C57BL/6 Mouse Model of SARS-CoV-2 Infection Recapitulates Age- and Sex-Based Differences in Human COVID-19 Disease and Recovery. Vaccines (Basel). 2022 Dec 25;11(1):47. doi: 10.3390/vaccines11010047. PMID: 36679892; PMCID: PMC9860616.
