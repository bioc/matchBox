\name{computeCat}

\alias{computeCat}

\title{Computing overlap proportions among ordered vectors}

\description{
  \code{computeCat} computes the overlap proportions between
  pairs of ordered vectors of identifiers.
  The input to this function is a data.frame containing non-redundant
  identifiers and a number of ranking statistics organized by columns.
  This function enables comparing all possible pair combinations,
  or selecting one column as the reference ranking for the remaining.
  The output of this function can be used as the input to
  \code{plotCat}, which creates correspondence at the top
    curves, as used in  Irizarry et al, Nat Methods (2005), for
    comparing differential gene expression across platforms and labs.
}
  
\usage{
computeCat(data, size=nrow(data), idCol=1, ref,
			method = c("equalRank", "equalStat"),
			decreasing = TRUE)
}

\arguments{
  \item{data}{A data.frame produced by \code{mergeData},
    containing a column of unique identifiers and at least two columns
    of ranking statistics (e.g. t-statistics, fold-change, Cox coefficients)}
  \item{size}{numeric. The number of top ranking statistics
    to be considered in the computation of the overlap proportions.
    If omitted all rows in \code{data} will be considered.
    If \code{size} is large computation time may be long.}
  \item{idCol}{numeric or character. The index (by default equal to one),
    or the name of the column containing the common identifiers
    (e.g. ENTREZID, SYMBOLS, \dots). }
  \item{ref}{character. The column name corresponding to
    the ranking statistics to be used as the reference in all pairs of
    comparisons.}
  \item{method}{character. The method used to compute the overlap
    proportion between two ordered vectors of identifiers: either "equalRank" 
    or "equalStat". The first method computed the overlap based
    on equal ranks, whereas the latter uses equal statistics. }
  \item{decreasing}{logical. This argument defines whether
    decreasing or increasing ordering should be used }
}

\details{
  \code{computeCat} computes overlapping proportions
  between pairs of ordered vectors of identifiers.
  This function first finds all possible pairs of vector combinations,
  then it computes the corresponding overlapping
  proportions. If a column is selected as the reference,
  using the argument \code{ref}, only the combinations
  involving this column will be returned.
  
  Briefly, for each CAT curve two vectors of identifiers
  are first ordered by the ranking statistics of choice,
  then the overlap between the two vectors is computed
  by considering more and more identifiers (vector size). 
  
  This function enables to compute overlapping proportions
  using two distinct methods: "equalRank" or "equalStat".
  With "equalRank" the overlap is obtained between vectors
  of the same size using equal ranks, which in turn can
  potentially correspond to ranking statistics of different
  magnitude (e.g. the vectors are of the same
  size, but might have different ranking statistics).
  With "equalStat" the overlap is obtained between vectors
  defined by using equal ranking statistics, which can
  potentially correspond to different rank, and hence to
  vectors of different size (e.g. the vectors are of different
  size, but have similar ranking statistics).
}

\value{
  A list of lists in which each element correspond to a
  CAT curve. If a specific reference column is provided
  through the \code{ref} argument, the number of
  list elements is equal to the number of combinations
  involving the reference group, otherwise  all  possible
  combinations are returned.
  When the "equalRank" method is used each list element
  contains only the overlapping proportion, while when
  the "equalStat" method is used the number of genes
  with equal statistics is stored along with the overlapping
  proportion.
  This output is used to produce CAT curves,
  using the \code{plotCat} function, as described
  in Irizarry et al, Nat Methods (2005).
}


\note{
  Given the combinatorial nature of the computation,
  a long computational time can be necessary if the input
  \code{data} contains many columns and many rows
  (number of features).
  In such a case consider limiting the number of rows
  used using the \code{size} argument.
}


\seealso{
  See \code{\link{mergeData}} and \code{\link{plotCat}}.
}


\references{
  Irizarry, R. A.; Warren, D.; Spencer, F.; Kim, I. F.; Biswal, S.;
  Frank, B. C.; Gabrielson, E.; Garcia, J. G. N.; Geoghegan, J.;
  Germino, G.; Griffin, C.; Hilmer, S. C.; Hoffman, E.;
  Jedlicka, A. E.; Kawasaki, E.; Martinez-Murillo, F.;
  Morsberger, L.; Lee, H.; Petersen, D.; Quackenbush, J.;
  Scott, A.; Wilson, M.; Yang, Y.; Ye, S. Q.
  and Yu, W. Multiple-laboratory comparison of microarray platforms.
  Nat Methods, 2005, 2, 345-350

  Ross, A. E.; Marchionni, L.; Vuica-Ross, M.; Cheadle, C.;
  Fan, J.; Berman, D. M.; and Schaeffer E. M.
  Gene Expression Pathways of High Grade Localized Prostate Cancer.
  Prostate 2011, 71, 1568-1578
  
  Benassi, B.; Flavin, R.; Marchionni, L.; Zanata, S.; Pan, Y.;
  Chowdhury, D.; Marani, M.; Strano, S.; Muti, P.; and Blandino, G.
  c-Myc is activated via USP2a-mediated modulation of microRNAs
  in prostate cancer. Cancer Discovery, 2012, March, 2, 236-247
  
}


\author{ Luigi Marchionni \email{marchion@jhu.edu}}


\examples{
###load data
data(matchBoxExpression)

###the column name for the identifiers
idCol <- "SYMBOL"

###the column name for the ranking statistics
byCol <- "t"

###use lapply to remove redundancy from all data.frames
###default method is "maxORmin"
newMatchBoxExpression <- lapply(matchBoxExpression, filterRedundant, idCol=idCol, byCol=byCol)

###select t-statistics and merge into a new data.frame using SYMBOL
mat <- mergeData(newMatchBoxExpression, idCol=idCol, byCol=byCol)

###Compute CAT for decreasing t-statistics: all genes
cpH2L <- computeCat(mat, idCol=1,decreasing=TRUE, method="equalRank")

###Compute CAT for increasing t-statistics:only the first 300 genes
cpL2H <- computeCat(mat, idCol=1, size=300, decreasing=FALSE, method="equalRank")

###Compute CAT for increasing t-statistics:only the first 300 genes
###use the second column as the reference
cpL2H.ref <- computeCat(mat, idCol=1, size=300, ref="dataSetA.t",
  decreasing=FALSE, method="equalRank")

}


\keyword{manip}
