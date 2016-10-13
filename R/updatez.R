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
