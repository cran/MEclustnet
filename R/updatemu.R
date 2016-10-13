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
