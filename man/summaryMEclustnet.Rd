\name{summaryMEclustnet}
\alias{summaryMEclustnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of MEclustnet object.
}
\description{
Summary of the output of the function MEclustnet which fits a mixture of experts latent position cluster model.
}
\usage{
summaryMEclustnet(fit, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{fit}{
An object storing the output of the function \code{\link{MEclustnet}}.
}
  \item{Y}{
The n x n binary adjacency matrix, with 0 down the diagonal, that was passed to \code{\link{MEclustnet}}.
}
}
\value{
A list with:
\item{AICM}{The value of the AICM criterion for the fitted model.}
\item{BICM}{The value of the BICM criterion for the fitted model.}
\item{BICMCMC}{The value of the BICMCMC criterion for the fitted model.}
\item{betamean}{The posterior mean vector of the regression coefficients for the link probabilities model.}
\item{betasd}{The standard deviation of the posterior distribution of beta.}
\item{taumean}{A matrix with G rows, detailing the posterior mean of the regression coefficients for the mixing proportions  model.}
\item{tausd}{The standard deviation of the posterior distribution of tau.}
\item{mumean}{A G x d matrix containing the posterior mean of the latent locations' mean.}
\item{meansd}{The standard deviation of the posterior distribution of mu.}
\item{sigma2mean}{A vector of length G containing the posterior mean of the latent locations' covariance.}
\item{sigma2sd}{The standard deviation of the posterior distribution of the latent locations' covariance.}
\item{Kmode}{A vector of length n detailing the posterior modal cluster membership for each node.}
\item{zmean}{An n x d matrix containing the posterior mean latent location for each node.}
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


\examples{
data(us.twitter.adjacency)
data(us.twitter.covariates)
link.vars = c(1)
mix.vars = c(1,5,7,8)
\dontrun{fit = MEclustnet(us.twitter.adjacency, us.twitter.covariates, 
link.vars, mix.vars, G=4, d=2, itermax = 500, burnin = 50, 
uphill = 1, thin=10)
plot(summ$zmean, col=summ$Kmode)
summ = summary(fit, us.twitter.adjacency)}
}

\keyword{ methods }