% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcloglikelihood.R
\name{calcloglikelihood}
\alias{calcloglikelihood}
\title{Calculate the log likelihood function of the data.}
\usage{
calcloglikelihood(pis, y)
}
\arguments{
\item{pis}{Vector of link probabilities.}

\item{y}{Vector version of the adjacency matrix, with the diagonal removed.}
}
\value{
The value of the log likelihood function.
}
\description{
This function calculates the log likelihood function of the data.
}
\references{
Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
}
\seealso{
\code{\link{MEclustnet}}
}
