#' Calculate the log likelihood function of the data.
#'
#'This function calculates the log likelihood function of the data.
#'
#' @param pis Vector of link probabilities.
#' @param y Vector version of the adjacency matrix, with the diagonal removed.
#'
#' @return The value of the log likelihood function.
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
calcloglikelihood <-
function(pis, y)
{
  sum((y*log(pis))+((1-y)*log(1-pis)))
}
