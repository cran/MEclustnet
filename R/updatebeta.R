#' Update the logistic regression parameters in the link probabilities model.
#'
#' @description The Metropolis-Hastings update step for the logistic regression parameters in the link probabilities model, using a surrogate proposal distribution.
#' @param beta Vector of regression coefficients in the link probabilities.
#' @param p Length of beta.
#' @param x.link Matrix, with \eqn{n^2 - n} rows and the same number of columns as covariates (including the intercept), giving the differences in covariates for all pairs of nodes.
#' @param delta Vector of Euclidean distances between locations in the latent space of all pairs of nodes.
#' @param y Vector version of the adjacency matrix, with the diagonal removed.
#' @param epsilon Mean of the multivariate normal prior on beta.
#' @param psi Covariance of the multivariate normal prior on beta.
#' @param psi.inv Inverse covariance of the multivariate normal prior on beta.
#' @param pis Vector of length \eqn{n^2 - n} providing the link probabilities between all pairs of nodes.
#' @param countbeta Counter for number of steps for which the proposed beta value was accepted.
#' @param rho Scaling factor to be used to adjust the acceptance rate.
#' @param n.tilde Length of the vector version of the adjacency matrix, with the diagonal removed i.e. \eqn{n^2 - n}.
#'
#' @details See appendix of the paper detailed below for details.
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
#' @return A list:
#' \describe{
#' \item{beta}{The returned version of the beta parameter vector.}
#' \item{countbeta}{The count of the number of acceptances of beta to that point in the MCMC chain.}
#' }
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom stats runif
updatebeta <-
function(beta, p, x.link, delta, y, epsilon, psi, psi.inv, pis, countbeta, rho, n.tilde)
{
  t.x.link.delta = matrix(c(x.link, delta), p+1, n.tilde, byrow=TRUE)
  beta.tilde = matrix(c(beta, 1), nrow=1)
  const = c(1+exp(beta.tilde%*%t.x.link.delta))     # Constant term in quad approx covariance
  propcov = rho*solve((t(x.link*sqrt(pis/const))%*%(x.link*sqrt(pis/const)))+psi.inv)   # quad approx
  propmu = colSums(as.matrix((x.link*(y-pis))))+(psi.inv%*%t(epsilon-beta))               # quad approx
  propmu = c(propcov%*%propmu)                                                   # quad approx

  betastar = beta+rmvnorm(1,propmu,propcov)
  pisstar = calcpis(betastar, x.link, delta, n.tilde)

  betastar.tilde = matrix(c(betastar, 1), nrow=1)
  conststar = c(1+exp(betastar.tilde%*%t.x.link.delta))                    # Constant term in quad approx covariance
  propcovnew = rho*solve((t(x.link*sqrt(pisstar/conststar))%*%(x.link*sqrt(pisstar/conststar)))+psi.inv)   # quad approx
  propmunew = colSums(as.matrix(x.link*(y-pisstar)))+(psi.inv%*%t(epsilon-betastar))                       # quad approx
  propmunew = c(propcov%*%propmu)                                                                   # quad approx


  logacceptbeta = calcloglikelihood(pisstar, y) + dmvnorm(betastar, epsilon, psi, log = TRUE) + dmvnorm(beta-betastar,propmunew,propcovnew,log = TRUE)- calcloglikelihood(pis, y) - dmvnorm(beta, epsilon, psi, log = TRUE) - dmvnorm(betastar-beta, propmu, propcov,log = TRUE)
  if(log(runif(1)) <=  min(logacceptbeta, 0))
  {
    beta = betastar; countbeta = countbeta+1
  }
  list(beta, countbeta)
}
