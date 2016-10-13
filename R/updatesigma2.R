updatesigma2 <-
function(G, alpha, m, d, sigma02, z, K, mu, sigma2)                  # Update variances in each cluster
{
  for(g in 1:G)
  {
    dof = alpha + (m[g]*d)
    if(sum(K == g)  ==  0){sc = sigma02}else{ sc = sigma02 + sum((t(z[(K == g),]) - mu[g,])^2)}
    sigma2[g] = (sc/rchisq(1, df = dof))
  }
  sigma2
}
