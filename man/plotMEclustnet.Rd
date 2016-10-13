\name{plotMEclustnet}
\alias{plotMEclustnet}

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
A vector detailing the columns of the covariate matrix which are included in the link probabilities regression model.
}
  \item{mix.vars}{
A vector detailing the columns of the covariate matrix which are included in the mixing proportions regression model.
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
clusters
}
}

\keyword{ hplot}
