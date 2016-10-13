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
