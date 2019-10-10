#' Label switching correction.
#'
#'This function corrects for the issue of label switching when fitting mixture models in a Bayesian setting.
#'
#' @param mu A G x d matrix of mean latent locations.
#' @param sigma2 A vector of length G containing the covariance of the latent locations within each cluster.
#' @param lambda An n x G matrix of mixing proportions.
#' @param tau A matrix of logistic regression coefficients, with G rows and number of columns equal to the number of covariates in the mixing proportions model plus 1, for the intercept.
#' @param K Vector of length n detailing the number of the cluster to which each node belongs.
#' @param G The number of clusters in the model being fitted.
#' @param d The dimension of the latent space.
#' @param perms A G! x G matrix of all possible permutations of 1:G (output by permutations(G), say).
#' @param muMAP A G x d matrix of \emph{maximum a posteriori} latent location means, obtained at the end of the uphill only section of the MCMC chain. Used as the template to correct for label switching.
#' @param iter Iteration number.
#' @param uphill Number of iterations for which uphill only steps in the MCMC chain should be run.
#' @param burnin Number of iterations of the MCMC chain which should not be included in \emph{a posteriori} summaries.
#' @param thin Thinning frequency of the MCMC chain to ensure independent samples.
#' @param s Number of columns in the reformatted covariates matrix for the mixing proportions model, output by \code{\link{formatting.covars}}.
#' @param x.mix The reformatted covariates matrix for the mixing proportions model, output by \code{\link{formatting.covars}.}
#'
#' @details The muMAP matrix is used as the reference to which each new estimate the cluster means is matched to correct for any label switching which may have occurred during sampling. A sum of squares function is employed as the loss function.
#' @return A list containing:list(mu, sigma2, lambda, tau, K)
#'\describe{
#'\item{mu}{The label-corrected matrix of cluster means.}
#'\item{sigma2}{The label-corrected vector of cluster covariances.}
#'\item{lambda}{The label-corrected matrix of mixing proportions.}
#'\item{tau}{The label-corrected matrix of logistic regression coefficients for the mixing proportions model.}
#'\item{K}{The label-corrected vector of length n detailing the number of the cluster to which each node belongs.}
#'}
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
labelswitch <-
function(mu, sigma2,lambda, tau, K, G, d, perms, muMAP, iter, uphill, burnin, thin, s, x.mix)
{
  smallloss = index = 0
  muperm = mu
  muold = mu
  sigma2old = sigma2
  tauold = tau
  Kold = K

  for(v in 1:factorial(G))
  {
    for(g in 1:G)
    {
      muperm[g,] = mu[perms[v, g],]
    }   #g
    loss = sum((muperm-muMAP)*(muperm-muMAP))
    if(v  ==  1)
    {
      smallloss = loss
      index = 1
    }
    if(loss < smallloss)
    {
      smallloss = loss
      index = v
    }
  }  #v

  # Relabel groups
  mu = matrix(muold[perms[index,],], G, d)
  sigma2 = sigma2old[perms[index,]]
  tau = matrix(tauold[perms[index,],], G, s)
  K = perms[index,Kold]

  # Ensure tau group 1 parameters are 0. Subsequently calculate lambda.
  tau = t(t(tau) - tau[1,])
  lambda = calclambda(tau, x.mix)

  list(mu, sigma2, lambda, tau, K)
}
