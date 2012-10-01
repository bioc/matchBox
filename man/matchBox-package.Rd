\name{matchBox-package}
\alias{matchBox-package}
\alias{matchBox}
\docType{package}

\title{
  A package to produce Correspondence-At-Top plots (CAT-plots) between ranked list of genes.
}

\description{
  The matchBox package allows to annotate and compare
  ranked vectors (e.g. by differential expression) of genomic features
  (e.g. genes, or probe sets), using CAT curves.
  A CAT curve displays the overlap proportion between two ranked
  vectors of identifiers against the number of considered features.
  This techiques was used for comparing differential gene expression
  results obtained from different platforms in different laboratories
  (see Irizarry et al, Nat Methods (2005))
}

\section{matchBox package features}{
  \describe{
    {Enables to filter data.frames containing feature
      identifiers and ranking statistics;}
    
    {Enables to identify the common set of features across
      results from different genomic experiments;}
    
    {Enables to merge multiple data.frames based on
      the commons set of features;}
    
    {Computes the overlap proportion between any pair of ranked
      features;}
    
    {Creates plots of the proportion of overlap along
      with conidence intervals;}
}}

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

\keyword{ package }

