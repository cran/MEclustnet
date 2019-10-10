#' Title Compute mixing proportions
#'
#'Function to compute the each observation's mixing proportions which are modeled as a logistic function of their covariates.
#'
#' @param tau A matrix of logistic regression coefficients, with G rows and number of columns equal to the number of covariates in the mixing proportions model plus 1, for the intercept.
#' @param x.mix A matrix of covariates in the mixing proportions model (including dummy variables for any factor covariates), with a column of 1's appended at the front.
#'
#' @return An n x G matrix of mixing proportions.
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
calclambda <-
function(tau, x.mix)
{
  res = tau%*%t(x.mix)
  res = t(exp(res-apply(res,2,max)))
  sweep(res, 1, rowSums(res), "/")
}
