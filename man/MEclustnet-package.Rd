\name{MEclustnet-package}
\alias{MEclustnet-package}
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
Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.}

\keyword{ package }
\keyword{ cluster }


\examples{
# An example from the Gormley and Murphy (2010) paper, using the Lazega lawyers friendship network.
# Iteration etc. are set to low values for illustrative purposes.
#Longer run times are likely to be required to achieve sufficient mixing.
data(lawyers.adjacency.friends)
data(lawyers.covariates)
link.vars = c(1)
mix.vars = c(1,4,5)
\dontrun{fit = MEclustnet(lawyers.adjacency.friends, lawyers.covariates,
link.vars, mix.vars, G=2, d=2, itermax = 500, burnin = 50, uphill = 1, thin=10)
# Plot the trace plot of the mean of dimension 1 for each cluster.
matplot(t(fit$mustore[,1,]), type="l")
# Compute posterior summaries
summ = summaryMEclustnet(fit, lawyers.adjacency.friends)
plot(summ$zmean, col=summ$Kmode, xlab="Dimension 1", ylab="Dimension 2", pch=summ$Kmode)
# Plot the resulting latent space, with uncertainties
plotMEclustnet(fit, lawyers.adjacency.friends)}

# An example analysing a 2016 Twitter network of US politicians.
# Iteration etc. are set to low values for illustrative purposes.
#Longer run times are likely to be required to achieve sufficient mixing.
data(us.twitter.adjacency)
data(us.twitter.covariates)
link.vars = c(1)
mix.vars = c(1)
\dontrun{fit = MEclustnet(us.twitter.adjacency, us.twitter.covariates,
link.vars, mix.vars, G=4, d=2, itermax = 500, burnin = 50, uphill = 1, thin=10)
# Plot the trace plot of the mean of dimension 1 for each cluster.
matplot(t(fit$mustore[,1,]), type="l")
# Compute posterior summaries
summ = summaryMEclustnet(fit, us.twitter.adjacency)
plot(summ$zmean, col=summ$Kmode, xlab="Dimension 1", ylab="Dimension 2", pch=summ$Kmode)
# Plot the resulting latent space, with uncertainties
plotMEclustnet(fit, us.twitter.adjacency)
# Examine which politicians are in which clusters...
clusters = list()
for(g in 1:fit$G)
{
  clusters[[g]] = us.twitter.covariates[summ$Kmode==g,c("name", "party")]
}
clusters}
}
