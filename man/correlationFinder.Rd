\name{correlationFinder}
\alias{correlationFinder}
\title{Find correlations in yeast expression data}
\usage{
  correlationFinder(dataFile = system.file("extdata", "sub_combined_complete_dataset_526G_198E.txt", package = "yeastCorrelationFinder"))
}
\arguments{
  \item{dataFile}{A tab-delimited yeast expression data
  file.}
}
\value{
  A named list where the names are a gene name, and the
  values are a list of genes which are correlated.
}
\description{
  Find genes which have a correlation between their
  expression profiles is equal to 0.84.
}
\details{
  Allocco (2004) stated, "In S. cerevisiae, two genes have
  a 50\% chance of having a common transcription factor
  binder if the correlation between their expression
  profiles is equal to 0.84." This function finds such
  correlations.
}
\examples{
correlationFinder()
}

