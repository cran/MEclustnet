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
