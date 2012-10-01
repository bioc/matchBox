\name{matchBoxExpression}

\alias{matchBoxExpression}

\docType{data}

\title{Example data: ranking from three differential gene expression experiments}

\description{
  List of differentially expressed genes from three distinct experiments
  from which the identifiers and the ranking statistics to be used
  for computing overlap proportions will be retrieved.
  This type of object is the starting point of a CAT-plot analysis,
}

\usage{data(matchBoxExpression)}

\format{
  This object is a list of data.frames containing at least two
  common columns, one for the identifiers and one for
  the ranking statistics. The common columns must
  have the same column names.
  In the  provided example the following columns
  are present in each data.frame:
  
  \describe{
    {\code{SYMBOL}: Gene symbol column;}
    
    {\code{GENENAME}: Gene name column;}

    {\code{ENTREZID}: ENTREZ Gene identifier column;}

    {\code{logFC}: Log2 fold-change column;}

    {\code{AveExpr}: Average expression (A-value) column;}

    {\code{t}: moderated t-statistics column;}

    {\code{P.Value}: P-value column;}

    {\code{adj.P.Val}: adjusted P-value column;}

    {\code{B}: B-statitics (log-odds) column;}
    
    The statistics were computed using the \code{topTable}
    function from \code{limma}.
  }
}


\author{
Luigi Marchionni \email{marchion@jhu.edu}
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

\keyword{datasets}
