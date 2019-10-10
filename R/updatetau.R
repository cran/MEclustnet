#' Update the logistic regression parameters in the mixing proportions model.
#'
#' @description The Metropolis-Hastings update step for the logistic regression parameters in the mixing proportions model, using a surrogate proposal distribution.
#'
#' @param G The number of clusters being fitted.
#' @param x.mix A matrix of covariates in the mixing proportions model (including dummy variables for any factor covariates), with a column of 1's appended at the front.
#' @param lambda An n x G matrix of mixing proportions.
#' @param Sigmag Covariance matrix of the multivariate normal prior for tau.
#' @param Sigmag.inv The inverse of Sigmag.
#' @param K The cluster membership vector.
#' @param gammag Mean vector of the multivariate normal prior for tau.
#' @param tau A matrix of logistic regression coefficients, with G rows and number of columns equal to the number of covariates in the mixing proportions model plus 1, for the intercept.
#' @param counttau Counter for number of steps for which the proposed tau value was accepted.
#' @param rho Scaling factor to be used to adjust the acceptance rate.
#'
#' @return A list:
#' \describe{
#' \item{tau}{The returned version of the tau parameter vector.}
#' \item{lambda}{The returned version of the lambda matrix.}
#' \item{counttau}{The count of the number of acceptances of tau to that point in the MCMC chain.}
#' }
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom stats runif
updatetau <-
function(G, x.mix, lambda, Sigmag, Sigmag.inv, K, gammag, tau, counttau, rho)
{
  if(G > 1)
  {
    for(g in 2:G)
    {
      x.lam.lam = x.mix*sqrt(lambda[,g]*(1 - lambda[,g]))
      propcov = rho*solve(t(x.lam.lam)%*%x.lam.lam + Sigmag.inv) # quadratic approx
      propmu = colSums(x.mix*((K == g)-(lambda[,g])))+(Sigmag.inv%*%t(gammag-tau[g,]))         # quadratic approximation
      propmu = c(propcov%*%propmu)

      change = rmvnorm(1,propmu,propcov)
      taustar = tau
      taustar[g,] = tau[g,]+change
      lambdastar = calclambda(taustar, x.mix)

      x.lamstar.lamstar = x.mix*sqrt(lambdastar[,g]*(1 - lambdastar[,g]))
      propcovnew = rho*solve(t(x.lamstar.lamstar)%*%x.lamstar.lamstar + Sigmag.inv)  # Quad approx
      propmunew = colSums(x.mix*((K == g)-(lambdastar[,g])))+(Sigmag.inv%*%t(gammag-taustar[g,]))          #quadratic approximation
      propmunew = c(propcovnew%*%propmunew)

      num = sum(((K == g)*log(lambdastar[,g]))+((1-(K == g))*log(lambdastar[,1]))) + dmvnorm(taustar[g,], gammag, Sigmag, log = TRUE) + dmvnorm(tau[g,]-taustar[g,],propmunew,propcovnew,log = TRUE)
      denom = sum(((K == g)*log(lambda[,g]))+((1-(K == g))*log(lambda[,1]))) + dmvnorm(tau[g,], gammag, Sigmag, log = TRUE) + dmvnorm(taustar[g,]-tau[g,],propmu,propcov,log = TRUE)
      logaccepttau = num - denom
      if (is.na(logaccepttau)){logaccepttau = 0}
      if(log(runif(1)) <=  min(logaccepttau, 0))
      {
        tau = taustar
        lambda = lambdastar
        counttau = counttau+1
      }
    } #g
  }
  list(tau, lambda, counttau)
}
