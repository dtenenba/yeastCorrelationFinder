## Run this before installing package, in the directory above the package:
## R --vanilla -e "library(roxygen2); roxygenize('yeastCorrelationFinder')"

##' Find correlations in yeast expression data
##'
##' Find genes which have a correlation between
##' their expression profiles is equal to 0.84.
##' @param dataFile A tab-delimited yeast expression data file.
##' @details Allocco (2004) stated, "In S. cerevisiae, 
##' two genes have a 50\% chance of having a common transcription
##' factor binder if the correlation between their expression
##' profiles is equal to 0.84." This function finds such correlations.
##' @return A named list where the names are a gene name,
##' and the values are a list of genes which are correlated.
##' @examples
##' correlationFinder()
##' @export
##' @importFrom gdata trim

correlationFinder <- function(dataFile = system.file("extdata",
    "sub_combined_complete_dataset_526G_198E.txt", 
    package="yeastCorrelationFinder"))
{
  tbl <- read.table (dataFile, sep='\t', header=T, quote='',
    comment.char='', fill=T, stringsAsFactors=FALSE)
  rownames (tbl) <- tbl$X
  exclude.these.columns <- 
  which (sapply (1:ncol (tbl),
    function (col) class (tbl [,col])) != 'numeric')
  if (length (exclude.these.columns) > 0)
    tbl <- tbl [, -exclude.these.columns]
  mtx.cor  <- cor (t (as.matrix (tbl)), use='pairwise.complete.obs')
  mtx.cor <- upper.tri (mtx.cor) * mtx.cor

  max = nrow (mtx.cor)
  correlated.genes <- list ()
  cor.threshold <- 0.85
  ret <- list()
  for (r in 1:max) {
    zz = as.integer (which (mtx.cor [r,] > cor.threshold))
    if (length (zz) > 0) {
      gene.a = rownames (mtx.cor) [r]
      genes.b = rownames (mtx.cor) [zz]
      correlated.genes [[gene.a]] <- genes.b
      ret[[ rownames(mtx.cor)[r] ]] <-
        rownames(mtx.cor)[zz]
  }
}
  ret
}
