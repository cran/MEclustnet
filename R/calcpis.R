#' Calculate link probabilities.
#'
#' Function calculates link probabilities between nodes.
#'
#' @param beta Vector of regression coefficients in the link probabilities.
#' @param x.link Matrix, with \eqn{n^2 - n} rows and the same number of columns as covariates (including the intercept), giving the differences in covariates for all pairs of nodes.
#' @param delta Vector of Euclidean distances between locations in the latent space of all pairs of nodes.
#' @param n.tilde Length of the vector version of the adjacency matrix, with the diagonal removed i.e. \eqn{n^2 - n}.
#'
#' @return A vector of length \eqn{n^2 - n} providing the link probabilities between all pairs of nodes.
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
calcpis <-
function(beta, x.link, delta, n.tilde)
{
  beta.tilde = matrix(c(beta,1), nrow=1)
  x.link.tilde = matrix(c(x.link, delta),ncol=n.tilde, byrow=T)
  res = exp(beta.tilde%*%x.link.tilde)
  c(res/(1 + res))
}
