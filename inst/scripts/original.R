## This is the original ad hoc script.
list.to.string = function(list.of.items, sep = " ")
{
    s = list.of.items[1]
    if (length(list.of.items) == 1)
        return(gdata::trim(s))
    for (item in list.of.items[2:length(list.of.items)])
    	s = paste(s, item, sep = sep)
    return(gdata::trim(s))
}

printf <- function(msg, ...) print(noquote(sprintf(msg, ...)))

f <- "sub_combined_complete_dataset_526G_198E.txt"
tbl <- read.table (f, sep='\t', header=T, quote='', comment.char='', fill=T, stringsAsFactors=FALSE)
rownames (tbl) <- tbl$X
exclude.these.columns <- 
which (sapply (1:ncol (tbl), function(col) class (tbl [,col])) != 'numeric')
if (length (exclude.these.columns) > 0)
  tbl <- tbl [, -exclude.these.columns]
mtx.cor  <- cor (t (as.matrix (tbl)), use='pairwise.complete.obs')
mtx.cor <- upper.tri (mtx.cor) * mtx.cor

max = nrow (mtx.cor)
correlated.genes <- list ()
cor.threshold <- 0.85
for (r in 1:max) {
  zz = as.integer (which (mtx.cor [r,] > cor.threshold))
  if (length (zz) > 0) {
    gene.a = rownames (mtx.cor) [r]
    genes.b = rownames (mtx.cor) [zz]
    correlated.genes [[gene.a]] <- genes.b
    printf ('%s: %s', rownames (mtx.cor)[r], list.to.string (rownames (mtx.cor)[zz]))}
}
