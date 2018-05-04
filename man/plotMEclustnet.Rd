\name{plotMEclustnet}
\alias{plotMEclustnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot latent position network.
}
\description{
Function to plot the resulting fitted network, using first two dimensions only.
}
\usage{
plotMEclustnet(fit, Y, link.vars, mix.vars)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
An object storing the output of the function \code{\link{MEclustnet}}.
}
  \item{Y}{
The n x n binary adjacency matrix, with 0 down the diagonal, that was passed to \code{\link{MEclustnet}}.
}
\item{link.vars}{
A vector of the column numbers of the data frame \code{covars} to be included in link probability model. If none are to be included, this argument should be 1.
}
  \item{mix.vars}{
A vector of the column numbers of the data frame \code{covars} to be included in mixing proportions model. If none are to be included, argument should be 1.
}
}
\details{
This function will plot the posterior mean latent location for each node in the network. The colour of each node reflects the posterior modal cluster membership, and the ellipses are 50\% posterior sets illustrating the uncertainty in the latent locations. The grey lines illustrate the observed links between the nodes.
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
summ = summary(fit, us.twitter.adjacency)
plot(summ$zmean, col=summ$Kmode)}
}

\keyword{ hplot}