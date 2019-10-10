#' Update step for the latent locations.
#'
#' @description A Metropolis-Hastings update step for the latent locations.
#'
#' @param n The number of nodes.
#' @param z The n x d matrix of latent locations.
#' @param x.link Matrix, with \eqn{n^2 - n} rows and the same number of columns as covariates (including the intercept), giving the differences in covariates for all pairs of nodes.
#' @param delta Vector of Euclidean distances between locations in the latent space of all pairs of nodes.
#' @param beta Vector of regression coefficients in the link probabilities.
#' @param y Vector version of the adjacency matrix, with the diagonal removed.
#' @param mu The G x d matrix of cluster means.
#' @param K The cluster membership vector
#' @param sigma2 The covariance of each cluster.
#' @param Id A d dimensional identity matrix.
#' @param pis A vector of length \eqn{n^2 - n} providing the link probabilities between all pairs of nodes.
#' @param iter Iteration number.
#' @param uphill Number of iterations for which uphill only steps in the MCMC chain should be run.
#' @param countz Counter for number of steps for which the proposed z value was accepted.
#' @param delete Index of the terms to be deleted in order to delete the diagonal terms from the vector version of the adjacency matrix.
#' @param d The dimension of the latent space.
#' @param n.tilde Length of the vector version of the adjacency matrix, with the diagonal removed i.e. \eqn{n^2 - n}.
#'
#' @return A list:
#' \describe{
#' \item{z}{The returned matrix of latent locations.}
#' \item{delta}{Vector of Euclidean distances between locations in the latent space of all pairs of nodes.}
#' \item{pis}{A vector of length \eqn{n^2 - n} providing the link probabilities between all pairs of nodes.}
#' \item{countz}{Counter for z acceptance rate.}
#' }
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
#' @importFrom mvtnorm dmvnorm
#' @importFrom MASS mvrnorm
#' @importFrom stats runif
updatez <-
function(n, z, x.link, delta, beta, y, mu, K, sigma2, Id, pis, iter,uphill, countz, delete, d, n.tilde)
{
  sigmaz2 = 1
  ind = sample(1:n, n, replace=F)
  for(j in 1:n)
  {
    i = ind[j]
    zstar = z
    zstar[i,] = mvrnorm(1, z[i,], sigmaz2*Id)
    deltastar = c(-as.matrix(dist(zstar)))[-delete]
    pisstar = calcpis(beta, x.link, deltastar, n.tilde)

    logacceptz = calcloglikelihood(pisstar, y) + dmvnorm(zstar[i,], mu[K[i],], sigma2[K[i]]*Id, log = T) - calcloglikelihood(pis, y) - dmvnorm(z[i,], mu[K[i],], sigma2[K[i]]*Id, log = T)
    if(iter < uphill)     # For first few runs force MH to accept uphill steps only to get MAP
    {
      if(logacceptz >=  0)
      {
        z = zstar
        delta = deltastar
        pis = pisstar
      } #if
    } #if
    if(iter >=  uphill)
    {
      if(log(runif(1)) <=  min(logacceptz, 0))
      {
        z = zstar
        delta = deltastar
        pis = pisstar
        countz = countz + 1
      }
    } # if
  } #j
  list(z, delta, pis, countz)
}
