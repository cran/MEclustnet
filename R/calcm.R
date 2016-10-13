calcm <-
function(m, G, K)                                               # Update counts in each cluster
{
  for(g in 1:G)
  {
    m[g]<-sum(K==g)
  } #g
  m
}
