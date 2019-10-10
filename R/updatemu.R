#' Update the mean of each cluster.
#'
#' @description A Gibbs step to update the mean of each cluster.
#' @param G The number of clusters being fitted.
#' @param z The n x d matrix of latent locations.
#' @param K The cluster membership vector.
#' @param m Vector of length G containing the number of nodes in each cluster.
#' @param sigma2 The covariance of each cluster.
#' @param omega2 Covariance of the multivariate normal prior distribution on the means. Note this is a scalar value, as the prior covariance is diagonal.
#' @param Id A d x d identity matrix.
#' @param mu The G x d matrix of cluster means.
#' @param d The dimension of the latent space.
#'
#' @return The G x d matrix of cluster means.
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
#' @importFrom MASS mvrnorm
updatemu <-
function(G, z, K, m, sigma2, omega2, Id, mu, d)                 # Update mean of each cluster
{
  for(g in 1:G)
  {
    if(m[g]  ==  0){meanvec = rep(0,d)}
    if(m[g] !=  0){meanvec = colSums(matrix(z[(K == g),], m[g], d))/(m[g] + (sigma2[g]/omega2))}
    covar = (sigma2[g]/(m[g] + (sigma2[g]/omega2)))*Id
    mu[g,] = mvrnorm(1, meanvec, covar)
  } #g
  mu
}
