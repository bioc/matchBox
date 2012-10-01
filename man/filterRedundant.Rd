\name{filterRedundant}

\alias{filterRedundant}


\title{This functions removes redundant features from a data.frame}

\description{
  Prior computing proportion of overlap between ranked vector of features
  it is necessary to remove the redundant features.
  This can be accomplished using a number of methods implemeted
  in the \code{filterRedundant} function, as explained below.
}

\usage{
filterRedundant(object,
    method=c("maxORmin", "geoMean", "mean", "median","random"),
    idCol=1, byCol=2, absolute=TRUE, decreasing=TRUE, trim=0, ...)
}

\arguments{
  \item{object}{a data.frame from which redundant features
    (rows) must be removed.}
  \item{method}{character. The method used for removing redundancy.
    Currently available methods are: \code{maxORmin}, 
    \code{geoMean}, \code{random}, \code{mean}, \code{median},
    (see Details below).}
  \item{idCol}{character or numeric. Name or index of the column
    containing redundant identifiers (e.g. ENTREZID, SYMBOLS, \dots).}
  \item{byCol}{character or numeric. Name or index of the column
    containing the ranking statistics (used only with \code{maxORmin}
    method).}
  \item{absolute}{logical. Indicates whether the absolute statistics,
    as defined by \code{byCol}, should be used when reordering
    (used only with \code{maxORmin} method).}
  \item{decreasing}{logical. Indicates whether reodering should be
    decreasing or not (used only with \code{maxORmin} method).}
  \item{trim}{numeric. Indicates whether a trimmed mean should
    be computed (used only with \code{mean} method).}
  \item{\dots}{further arguments to be passed (not currently implemented).}
}
 
\details{
  The \code{maxORmin} method removes
  redundant features by selecting the rows
  that correspond to the maximum or minimum
  value of a selected statistics.
  With this approach
  redundant features are first
  ranked in increasing or decreasing order,
  as defined by the \code{decreasing} argument,
  using the ranking statistics defined by \code{byCol},
  either in their original or absolute scale,
  as defined by \code{absolute} argument.
  Subsequently data.frame rows corresponding to redundant
  identifiers are removed, after these have been identified in
  the column defined by the \code{idCol},
  using the \code{duplicated} function.
  
  The \code{mean}, \code{median}, \code{geoMean},
  and \code{random} methods provide alternative ways
  for summarizing numerical values corresponding to
  redundant features, as defined by the \code{idCol}
  argument:
  \code{mean} takes the average,
  \code{median} the median,
  \code{geoMean} the geometric mean,
  \code{random} select a random value.
  
%%%   The \code{maxVar} method computes the
%%%   variance by rows, and uses it to rank the
%%%   features in decreasing order.
%%%   Redundancy is removed by keeping the feature
%%%   showing the greatest variance
%%%   (for details see Note section below).
  
}

\value{
  A data.frame with fewer rows with respect to the input one,
  unique by the identifier specified by the \code{idCol} argument.
}


\note{
  \code{filterRedundant} is a utility function providing various
  methods to remove redundant rows from a data.frame.
  The choice of the method depends on the nature of the values,
  and the final goal.
  %%% For instance if the data.frames contain
  %%% different statistics from a differential gene expression analysis
  %%% (like t-statistics, p-values, fold-changes) it would make no sense
  %%% to use the \code{maxVar} approach, which uses row variances.
  Therefore caution should be used when taking the mean
  or the median across few values, or passing the arguments
  with the \code{minORmax} method (for instance it would
  make no sense at all to use a decreasing ordering if the ranking
  statistics is a p-value).  
}


\seealso{
  See \code{\link[base]{duplicated}}.
}


\author{Luig Marchionni <marchion@jhu.edu>}

\examples{
###load data
data(matchBoxExpression)

###check whether there are redundant identifiers
sapply(matchBoxExpression,nrow)

###the column name for the identifiers
idCol <- "SYMBOL"

###the column name for the ranking statistics
byCol <- "t"

###use lapply to remove redundancy from all data.frames
###default method is "maxORmin"
newMatchBoxExpression <- lapply(matchBoxExpression, filterRedundant, idCol=idCol, byCol=byCol)

###recheck number of rows
sapply(newMatchBoxExpression, nrow)

}

\keyword{ manip }
