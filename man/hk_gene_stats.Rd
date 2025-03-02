% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hkgeneselection.r
\name{hk_gene_stats}
\alias{hk_gene_stats}
\title{house keeping gene analysis}
\usage{
hk_gene_stats(
  data,
  group.by = "sample",
  min_num_hk_genes = 5,
  manually_selected_hk_genes = FALSE,
  output_dir = getwd()
)
}
\arguments{
\item{data}{counts compiled from rcc files in load_nCounter_files step}

\item{group.by}{what group counts of each house keeping gene, individual samples are recomme5ded}

\item{min_num_hk_genes}{minimum number of house keeping genes to keep for analysis}

\item{manually_selected_hk_genes}{list of house keeping genes manually selected}

\item{output_dir}{directory for results (default is working directory)}
}
\description{
This assesses which genes should be used for ratio normalization
}
\examples{
data <- hk_gene_stats(counts)


}
\keyword{genes}
\keyword{house}
\keyword{keeping}
