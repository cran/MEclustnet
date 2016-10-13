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
