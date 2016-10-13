\name{updatetau}
\alias{updatetau}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update the logistic regression parameters in the mixing proportions model.
}
\description{
The Metropolis-Hastings update step for the logistic regression parameters in the mixing proportions model, using a surrogate proposal distribution.
}
\usage{
updatetau(G, x.mix, lambda, Sigmag, Sigmag.inv, K, gammag, tau, counttau, rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
The number of clusters being fitted.
}
  \item{x.mix}{
A matrix of covariates in the mixing proportions model (including dummy variables for any factor covariates), with a column of 1's appended at the front.
}
  \item{lambda}{
An n x G matrix of mixing proportions.
}
  \item{Sigmag}{
Covariance matrix of the multivariate normal prior for tau.
}
  \item{Sigmag.inv}{
The inverse of Sigmag.
}
  \item{K}{
The cluster membership vector.
}
  \item{gammag}{
Mean vector of the multivariate normal prior for tau.
}
  \item{tau}{
A matrix of logistic regression coefficients, with G rows and number of columns equal to the number of covariates in the mixing proportions model plus 1, for the intercept.
}
  \item{counttau}{
Counter for number of steps for which the proposed tau value was accepted.
}
  \item{rho}{
Scaling factor to be used to adjust the acceptance rate.
}
}

\value{
A list:
\item{tau}{The returned version of the tau parameter vector.}
\item{lambda}{The returned version of the lambda matrix.}
\item{counttau}{The count of the number of acceptances of tau to that point in the MCMC chain.}
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
