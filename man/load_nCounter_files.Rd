% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_nCounter_data.r
\name{load_nCounter_files}
\alias{load_nCounter_files}
\title{Load nCounter files}
\usage{
load_nCounter_files(
  pathtoRCC = getwd(),
  meta.data = "",
  save.fig = TRUE,
  output_dir = getwd()
)
}
\arguments{
\item{pathtoRCC}{matrix raw or normalized (default is current working directory)}

\item{meta.data}{file with experimental information (meta data csv file)}
}
\description{
This loads count and target file info (internal function only)
}
\examples{
data <- load_nCounter_files(pathtoRCC="./", meta.data="targetfile.csv")

}
\keyword{Ncounter}
\keyword{RCC}
\keyword{files}
