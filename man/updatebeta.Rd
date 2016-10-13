\name{updatebeta}
\alias{updatebeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update the logistic regression parameters in the link probabilities model.
}
\description{
The Metropolis-Hastings update step for the logistic regression parameters in the link probabilities model, using a surrogate proposal distribution.
}
\usage{
updatebeta(beta, p, x.link, delta, y, epsilon, psi, psi.inv, pis, countbeta, rho, n.tilde)
}
\arguments{
  \item{beta}{
Vector of regression coefficients in the link probabilities.
}
  \item{p}{
Length of beta.
}
  \item{x.link}{
Matrix, with \eqn{n^2 - n} rows and the same number of columns as covariates (including the intercept), giving the differences in covariates for all pairs of nodes.

}
  \item{delta}{
Vector of Euclidean distances between locations in the latent space of all pairs of nodes.
}
  \item{y}{
Vector version of the adjacency matrix, with the diagonal removed.
}
  \item{epsilon}{
Mean of the multivariate normal prior on beta.
}
  \item{psi}{
Covariance of the multivariate normal prior on beta.
}
  \item{psi.inv}{
Inverse covariance of the multivariate normal prior on beta.
}
  \item{pis}{
Vector of length \eqn{n^2 - n} providing the link probabilities between all pairs of nodes.
}
  \item{countbeta}{
Counter for number of steps for which the proposed beta value was accepted.
}
  \item{rho}{
Scaling factor to be used to adjust the acceptance rate.
}
  \item{n.tilde}{
Length of the vector version of the adjacency matrix, with the diagonal removed i.e. \eqn{n^2 - n}.
}
}
\details{
See appendix of the paper detailed below for details.
}
\value{
A list:
\item{beta}{The returned version of the beta parameter vector.}
\item{countbeta}{The count of the number of acceptances of beta to that point in the MCMC chain.}
%% ...
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
