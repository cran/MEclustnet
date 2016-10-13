\name{updatez}
\alias{updatez}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update step for the latent locations.
}
\description{
A Metropolis-Hastings update step for the latent locations.
}
\usage{
updatez(n, z, x.link, delta, beta, y, mu, K, sigma2, Id, pis,
iter, uphill, countz, delete, d, n.tilde)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
The number of nodes.
}
  \item{z}{
The n x d matrix of latent locations.
}
  \item{x.link}{
Matrix, with \eqn{n^2 - n} rows and the same number of columns as covariates (including the intercept), giving the differences in covariates for all pairs of nodes.

}
  \item{delta}{
Vector of Euclidean distances between locations in the latent space of all pairs of nodes.
}
  \item{beta}{
Vector of regression coefficients in the link probabilities.
}
  \item{y}{
Vector version of the adjacency matrix, with the diagonal removed.
}
  \item{mu}{
The G x d matrix of cluster means.
}
  \item{K}{
The cluster membership vector.
}
  \item{sigma2}{
The covariance of each cluster.
}
  \item{Id}{
A d dimensional identity matrix.
}
  \item{pis}{
A vector of length \eqn{n^2 - n} providing the link probabilities between all pairs of nodes.
}
  \item{iter}{
Iteration number.
}
  \item{uphill}{
Number of iterations for which uphill only steps in the MCMC chain should be run.
}
  \item{countz}{
Counter for number of steps for which the proposed z value was accepted.
}
  \item{delete}{
Index of the terms to be deleted in order to delete the diagonal terms from the vector version of the adjacency matrix.
}
  \item{d}{
The dimension of the latent space.
}
  \item{n.tilde}{
Length of the vector version of the adjacency matrix, with the diagonal removed i.e. \eqn{n^2 - n}.
}
}

\value{
 list(z, delta, pis, countz)
A list:
\item{z}{The returned matrix of latent locations.}
\item{delta}{Vector of Euclidean distances between locations in the latent space of all pairs of nodes.}
\item{pis}{A vector of length \eqn{n^2 - n} providing the link probabilities between all pairs of nodes.}
\item{countz}{Counter for z acceptance rate.}
}
\references{
Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
}
\author{
Isobel Claire Gormley <claire.gormley@ucd.ie>
}

\seealso{
\code{\link{MEclustnet}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
