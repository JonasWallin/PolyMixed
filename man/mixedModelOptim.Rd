% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneModel.R
\name{mixedModelOptim}
\alias{mixedModelOptim}
\title{estimate the gene model}
\usage{
mixedModelOptim(
  GeneMixObj,
  find = NULL,
  dupl = NULL,
  estPar = TRUE,
  restricted = T
)
}
\arguments{
\item{find}{- (k x 1) known non zero genes}

\item{dupl}{- (c x 1) duplicated rows of X (found in GeneMixObj$LikObj$X)}
}
\description{
estimate the gene model
}
