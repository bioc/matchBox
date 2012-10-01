\name{calcHypPI}

\alias{calcHypPI}

\title{Probability intervals calculation for CAT curves using the hypergeometric distribution.}

\description{
  The \code{calcHypPI} function calculates probability intervals
   for a correspondence at the top (CAT) curve using the
   hypergeometric distribution. This function, based on the
   \code{qhyper} quantile function,  produces a probability
   intervals matrix to be passed as argument to \code{plotCat}
   in order to add probability intervals shades when plotting CAT curves.
}

\usage{
calcHypPI(data,  expectedProp = 0.1, prob = c(0.999999,0.999,0.99,0.95))
}

\arguments{
  \item{data}{The same data frame used to compute the CAT curves
    with the \code{computeCat} function. It contains a column
    of unique identifiers and at least two columns of ranking statistics.}
  \item{expectedProp}{A single numeric value between 0 and 1.
    This is the proportion of features expected to be corresponding
    at the top of the ranking. The "expectedProp" argument can be set to NULL
    if the number of features expected to be similarly ranked is unknown.}
  \item{prob}{A numeric vector specifying the probabiliy intervals
    for the CAT curves to be computed.}
}

\details{
  The \code{calcHypPI} uses \code{qhyper} quantile function
  to compute the proportions of common features between two ordered
  vectors for specified quantiles of the hypergeometric distribution.
  Such proportions are used to add probability intervals
  to CAT curves computed using ranks (see \code{computeCat}).
  The \code{prob} argument is used to specify the desired probability
  intervals to be computed. By default this numeric vector is equal to
  \code{c(0.999999, 0.999, 0.99, 0.95)}.
 
  To understand the way this function works we can use
  the analogy of repeated drawing of an increasing number
  of balls from an urn containing both white and black balls
  (see \code{qhyper}).
  According to this analogy the total number of balls in the urn
  corresponds to the total number of common features 
  between two ordered vectors that are being compared
  (e.g. all the genes in common between two genomic studies).
  
  The number of white balls corresponds to the top ranking
  features that are correctly ordered (successes),
  while the black balls represent the features that are
  not correctly ordered (failures).

  Finally, according to this analogy, comparing the first top
  10 features from each vector will correspond to a first draw
  of 10 balls from the urn, while comparing the top 20
  features to a draw of 20 balls, and so on until all balls
  are drawn at once.

  By default the \code{calcHypPI} function expects
  that the top 10\% of the features of the two vectors
  are similarly ordered. This expectation can be modified
  by the \code{expectedProp} argument. When 
  \code{expectedProp} is set equal to \code{NULL}
  the number of white balls in the urn
  (i.e. the top ranking features in the correct order)
  corresponds to the number of balls that are drawn
  at each attempt (i.e. the increasing size of top features
  from each vector that are being compared).
}

\value{
  It returns a numeric matrix containing the probability intervals
  for CAT curves based on equal ranks.
  The column names of this matrix specifies the quantiles
  of the hypergeometric distribution used to compute
  the intervals. The values represent the proportions of overlap
  associated with the defined quantiles.
  The resulting matrix object is used to add the probability
  intervals shades when plotting CAT curves by passing it
  to the \code{preComputedPI} argument of the
  \code{plotCat} function. 
}

\note{
  This function will take more and more time to run when more
  and more features are used. For this reason it is convenient
  to compute the probability intervals separately and store
  the probability intervals matrix for re-use when plotting
  the CAT curves.
}

\seealso{
  See \code{\link[stats]{qhyper}}, \code{\link{plotCat}},
  \code{\link{calcHypPI}} and \code{\link{computeCat}}.
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
  Prostate, 2011, 71, 1568-1578

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

### compute probability intervals with default values
confInt <- calcHypPI(data=mat)

###structure of confInt
str(confInt)

### compute probability intervals with "expectedProp" set to NULL
confInt2 <- calcHypPI(data=mat, expectedProp=NULL)

###structure of confInt
str(confInt2)
}

\keyword{ manip }


