#' Totals the number of observations in each cluster.
#'
#'Update the count of the number of observations in each cluster.
#'
#' @param m Vector of length G containing the number of nodes in each cluster.
#' @param G The number of clusters in the model being fitted.
#' @param K Vector of length n detailing the number of the cluster to which each node belongs.
#'
#' @return Vector of length G containing the number of nodes in each cluster.
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
calcm <-
function(m, G, K)                                               # Update counts in each cluster
{
  for(g in 1:G)
  {
    m[g]<-sum(K==g)
  } #g
  m
}
