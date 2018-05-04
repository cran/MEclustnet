\name{MEclustnet-package}
\alias{MEclustnet-package}
%\alias{MEclustnet}
\docType{package}
\title{
\packageTitle{MEclustnet}
}
\description{
\packageDescription{MEclustnet}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{MEclustnet}
\packageIndices{MEclustnet}
The main function of interest is MEclustnet which will fit a mixture of experts latent position cluster model to a binary network.
}
\author{
\packageAuthor{MEclustnet}

Maintainer: \packageMaintainer{MEclustnet}
}
\references{
Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
}

\keyword{ package }
\keyword{ cluster }


\examples{

data(us.twitter.adjacency)
data(us.twitter.covariates)
link.vars = c(1)
mix.vars = c(1,5,7,8)
\dontrun{fit = MEclustnet(us.twitter.adjacency, us.twitter.covariates, link.vars, 
mix.vars, G=4, d=2, itermax = 500, burnin = 50, uphill = 1, thin=10)
summ = summary(fit, us.twitter.adjacency)
plot(summ$zmean, col=summ$Kmode) }
}
