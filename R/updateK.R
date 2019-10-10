#' Update the cluster membership vector.
#'
#' @description A Gibbs update step for K, the cluster membership vector.
#' @param G The number of clusters being fitted.
#' @param K The cluster membership vector.
#' @param z The n x d matrix of latent locations.
#' @param mu The G x d matrix of cluster means.
#' @param sigma2 The G vector of cluster covariances.
#' @param Id An identity matrix of dimension d.
#' @param lambda The n x G matrix of mixing proportions.
#'
#' @return The cluster membership vector.
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats rmultinom
updateK <-
function(G, K, z, mu, sigma2, Id, lambda)                     # Gibbs multinomial update step for K
{
  dens = rep(0,G)
  for(i in 1:length(K))
  {
    for(g in 1:G)
    {
      dens[g] = dmvnorm(z[i,], mu[g,], sigma2[g]*Id)
    } #g
    if(sum(is.na(dens))!= 0){dens[c(1:G)[is.na(dens)]] = 0}            # Check in case of empty group
    K[i] = c(1:G)[rmultinom(1,1,lambda[i,]*dens) == T]
  } #i
  K
}
