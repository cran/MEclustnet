#' Summary of MEclustnet object.
#'
#' @description Summary of the output of the function MEclustnet which fits a mixture of experts latent position cluster model.
#' @param fit An object storing the output of the function \code{\link{MEclustnet}}.
#' @param Y The n x n binary adjacency matrix, with 0 down the diagonal, that was passed to \code{\link{MEclustnet}}.
#'
#' @return A list with:
#' \describe{
#' \item{AICM}{The value of the AICM criterion for the fitted model.}
#' \item{BICM}{The value of the BICM criterion for the fitted model.}
#' \item{BICMCMC}{The value of the BICMCMC criterion for the fitted model.}
#' \item{betamean}{The posterior mean vector of the regression coefficients for the link probabilities model.}
#' \item{betasd}{The standard deviation of the posterior distribution of beta.}
#' \item{taumean}{A matrix with G rows, detailing the posterior mean of the regression coefficients for the mixing proportions  model.}
#' \item{tausd}{The standard deviation of the posterior distribution of tau.}
#' \item{mumean}{A G x d matrix containing the posterior mean of the latent locations' mean.}
#' \item{meansd}{The standard deviation of the posterior distribution of mu.}
#' \item{sigma2mean}{A vector of length G containing the posterior mean of the latent locations' covariance.}
#' \item{sigma2sd}{The standard deviation of the posterior distribution of the latent locations' covariance.}
#' \item{Kmode}{A vector of length n detailing the posterior modal cluster membership for each node.}
#' \item{zmean}{An n x d matrix containing the posterior mean latent location for each node.}
#' }
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
#' @importFrom stats var sd
#' @examples #################################################################
#' # An example analysing a 2016 Twitter network of US politicians.
#' #################################################################
#' # Number of iterations etc. are set to low values for illustrative purposes.
#' # Longer run times are likely to be required to achieve sufficient mixing.
#'
#' library(latentnet)
#' data(us.twitter.adjacency)
#' data(us.twitter.covariates)
#'
#' link.vars = c(1)
#' mix.vars = c(1,5,7,8)
#'
#' \donttest{fit = MEclustnet(us.twitter.adjacency, us.twitter.covariates,
#'  link.vars, mix.vars, G=4, d=2, itermax = 500, burnin = 50, uphill = 1, thin=10)
#'
#' # Plot the trace plot of the mean of dimension 1 for each cluster.
#' matplot(t(fit$mustore[,1,]), type="l", xlab="Iteration", ylab="Parameter")
#'
#' # Compute posterior summaries
#' summ = summaryMEclustnet(fit, us.twitter.adjacency)
#'
#' plot(summ$zmean, col=summ$Kmode, xlab="Dimension 1", ylab="Dimension 2", pch=summ$Kmode,
#'      main = "Posterior mean latent location for each node.")
#'
#' # Plot the resulting latent space, with uncertainties
#' plotMEclustnet(fit, us.twitter.adjacency, link.vars, mix.vars)
#'
#' # Examine which politicians are in which clusters...
#' clusters = list()
#' for(g in 1:fit$G)
#' {
#'   clusters[[g]] = us.twitter.covariates[summ$Kmode==g,c("name", "party")]
#' }
#' clusters
#' }
#' @export
summaryMEclustnet <-
function(fit, Y)
{
  # Compute AICM
  N = nrow(Y)
  meanll = mean(fit$LLstore,na.rm=TRUE)
  varll = var(fit$LLstore,na.rm=TRUE)
  AICM = 2*(meanll - varll)
  BICM = meanll - varll*(log(N) - 1)
  nu = ncol(fit$betastore) + (fit$G*fit$d) + fit$G + (fit$G*dim(fit$taustore)[2])
  BICMCMC = 2*max(fit$LLstore) - nu*log(N)

  betamean = apply(fit$betastore,2,mean,na.rm=TRUE)
  betasd = apply(fit$betastore,2,sd,na.rm=TRUE)


  taumean = apply(fit$taustore, c(1,2), mean, na.rm=T)
  tausd = apply(fit$taustore, c(1,2), sd, na.rm=T)


  mumean = apply(fit$mustore, c(1,2), mean, na.rm=T)
  musd = apply(fit$mustore, c(1,2), sd, na.rm=T)

  sigma2mean = apply(fit$sigma2store, 2, mean, na.rm=T)
  sigma2sd = apply(fit$sigma2store, 2, sd, na.rm=T)


  Kmode = as.numeric(apply(fit$Kstore,2,function(v){names(sort(-table(v)))[1]}))
  zmean = apply(fit$zstore,c(1:fit$d),mean,na.rm=TRUE)

  list(AICM = AICM, BICM = BICM, BICMCMC = BICMCMC, betamean = betamean, betasd = betasd, taumean = taumean, tausd = tausd, mumean = mumean, musd = musd, sigma2mean = sigma2mean, sigma2sd = sigma2sd, Kmode = Kmode, zmean = zmean)
}
