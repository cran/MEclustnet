#' Update variances in each cluster.
#'
#' @description A Gibbs step to update variances in each cluster.
#' @param G The number of clusters being fitted.
#' @param alpha Degrees of freedom of the scaled inverse Chi squared prior distribution on the cluster variances.
#' @param m Vector of length G containing the number of nodes in each cluster.
#' @param d Dimension of the latent space.
#' @param sigma02 Scaled factor of the scaled inverse Chi squared prior distribution on the cluster variances.
#' @param z The n x d matrix of latent locations.
#' @param K The cluster membership vector.
#' @param mu The G x d matrix of cluster means.
#' @param sigma2 The G vector of cluster variances.
#'
#' @return The G vector of cluster variances.
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
#' @importFrom stats rchisq
updatesigma2 <-
function(G, alpha, m, d, sigma02, z, K, mu, sigma2)                  # Update variances in each cluster
{
  for(g in 1:G)
  {
    dof = alpha + (m[g]*d)
    if(sum(K == g)  ==  0){sc = sigma02}else{ sc = sigma02 + sum((t(z[(K == g),]) - mu[g,])^2)}
    sigma2[g] = (sc/rchisq(1, df = dof))
  }
  sigma2
}
