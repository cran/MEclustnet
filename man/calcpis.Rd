% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcpis.R
\name{calcpis}
\alias{calcpis}
\title{Calculate link probabilities.}
\usage{
calcpis(beta, x.link, delta, n.tilde)
}
\arguments{
\item{beta}{Vector of regression coefficients in the link probabilities.}

\item{x.link}{Matrix, with \eqn{n^2 - n} rows and the same number of columns as covariates (including the intercept), giving the differences in covariates for all pairs of nodes.}

\item{delta}{Vector of Euclidean distances between locations in the latent space of all pairs of nodes.}

\item{n.tilde}{Length of the vector version of the adjacency matrix, with the diagonal removed i.e. \eqn{n^2 - n}.}
}
\value{
A vector of length \eqn{n^2 - n} providing the link probabilities between all pairs of nodes.
}
\description{
Function calculates link probabilities between nodes.
}
\references{
Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
}
\seealso{
\code{\link{MEclustnet}}
}
