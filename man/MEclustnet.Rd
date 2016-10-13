\name{MEclustnet}
\alias{MEclustnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a mixture of experts latent position cluster model to network data.
}
\description{
MEclustnet will fit a mixture of experts latent position cluster model to a binary network.
}
\usage{
MEclustnet(Y, covars, link.vars = c(1:ncol(covars)), mix.vars = c(1:ncol(covars)),
G = 2, d = 2, itermax = 10000, uphill = 100, burnin = 1000, thin = 10,
rho.input = 1, verbose = TRUE, ...)
}
\arguments{
  \item{Y}{
An n x n binary matrix of links between n nodes, with 0 on the diagonal and 1 indicating a link.
}
  \item{covars}{
An n x p data frame of node specific covariates. Categorical variables should be factors. First column should be a column of 1s, and should always be passed in.
}
  \item{link.vars}{
A vector of the column numbers of the data frame \code{covars} to be included in link probability model. If none are to be included, this argument should be 1.

}
  \item{mix.vars}{
A vector of the column numbers of the data frame \code{covars} to be included in mixing proportions model. If none are to be included, argument should be 1.

}
  \item{G}{
The number of clusters in the model to be fitted.
}
  \item{d}{
The dimension of the latent space.
}
  \item{itermax}{
Maximum number of iterations in the MCMC chain.
}
  \item{uphill}{
Number of iterations for which uphill only steps in the MCMC chain should be run to find \emph{maximum a posteriori} estimates.
}
  \item{burnin}{
Number of burnin  iterations in the MCMC chain.
}
  \item{thin}{
The degree of thinning to be applied to the MCMC chain.
}
  \item{rho.input}{
Scaling factor to achieve desireable acceptance rates in Metropolis-Hastings steps.
}
  \item{verbose}{
Print progress updates to screen? Recommended as the models are slow to run.
}
  \item{\dots}{
Additional arguments.
}
}
\details{
This function fits the mixture of experts latent position cluster model to a binary network via a Metropolis-within-Gibbs sampler. Covariates can influence either the link probabilities between nodes and/or the cluster memberships of nodes.
}
\value{
An object of class \code{MEclustnet}, which is a list containing:
\item{zstore}{An n x d x store.dim array of sampled latent location matrices, where store.dim is the number of post burnin thinned iterations.}
\item{betastore}{A store.dim x p matrix of sampled beta vectors, the logistic regression parameters of the link probabilities model.}
\item{Kstore}{A store.dim x n matrix of sampled cluster membership vectors.}
\item{mustore}{A G x d x store.dim array of sampled cluster mean latent location matrices.}
\item{sigma2store}{A store.dim x G matrix of sampled cluster variances.}
\item{lambdastore}{An n x G x store.dim array of sampled mixing proportion matrices.}
\item{taustore}{A G x s x store.dim array of sampled tau vectors, the logistic regression parameters of the mixing proportions model, where s is the length of tau.}
\item{LLstore}{A vector of length store.dim storing the loglikelihood from each stored iteration.}
\item{G}{The number of clusters fitted}
\item{d}{The dimension of the latent space}
\item{countbeta}{Count of accepted beta values}
\item{counttau}{Count of accepted tau values}
}
\references{
Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
}
\author{
Isobel Claire Gormley <claire.gormley@ucd.ie>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{summaryMEclustnet}}, \code{\link{plotMEclustnet}}
}
\examples{
# An example from the Gormley and Murphy (2010) paper, using the Lazega lawyers friendship network.
# Iteration etc. are set to low values for illustrative purposes.
# Longer run times are likely to be required to achieve sufficient mixing.
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
# Longer run times are likely to be required to achieve sufficient mixing.
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

\keyword{ methods }
\keyword{ cluster }
