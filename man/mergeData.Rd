\name{mergeData}

\alias{mergeData}

\title{Merging data.frames based on common identifiers}

\description{
  This utility function is used for merging specific columns
  from a set of distinct data.frames based on a specific set
  of identifiers. For instance this utility function can be
  used to retrieve from multiple data.frames the ranking
  statistics and the identifiers that will be used for computing
  the correspondence  at the top curves.
}

\usage{
mergeData(listOfDataFrames, idCol=1, byCol=2)
}

\arguments{
  \item{listOfDataFrames}{list. This object is a list
    of distinct data.frames to be merged based on
    common identifiers. The data.frames to be merged
    must contain at least  two common columns,
    one for the identifiers (as specified by \code{idCol}),
    and one for the ranking statistics (as specified by
    \code{byCol}).
    Redundant features are not allowed, and should
    be previously removed using \code{filterRedundant}. }
  \item{idCol}{character or numeris. Name or index of the column
    containing the common identifiers (e.g. ENTREZID, SYMBOLS, \dots).}
  \item{byCol}{character or numeric . Name of index the column
    containing the ranking statistics.}
}

\details{
  This function first identifies the common set of features
  across all the data.frames contained in the \code{listOfDataFrames}
  object. Subsequently, for this common set of features,
  it returns a single data.frame containing the ranking statistics
  values of choice collected from each data.frame.
}


\value{
  A data.frame containing the identifiers and
  the ranking statistics common to all data.frames
  in \code{listOfDataFrames} to be used for computing
  the correspondence at the top
  (see Irizarry et al, Nat Methods (2005))
}


\seealso{
  See \code{\link{filterRedundant}}.
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
mat <- mergeData(listOfDataFrames = newMatchBoxExpression, idCol = idCol,
byCol = byCol)

###structure of mat
str(mat)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }


