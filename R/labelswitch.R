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
