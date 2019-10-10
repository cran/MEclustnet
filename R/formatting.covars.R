#' Reformat matrix of covariates.
#'
#' This function reformats the matrix of input covariates into the required format for the link probabilities and for the mixing proportions.
#'
#' @param covars The n x p data frame of node specific covariates passed in to the overall \code{\link{MEclustnet}} function. The first column should be a column of 1's and categorical variables should be factors.
#' @param link.vars A vector detailing the column numbers of the matrix covars that should be included in the link probabilities model.
#' @param mix.vars A vector detailing the column numbers of the matrix covars that should be included in the mixing proportions probabilities model.
#' @param n The number of nodes in the network.
#'
#' @details For the link regression model, the difference in the link.vars covariates, for all pairs of nodes is calculated.
#' For the mixing proportions model, the required representation of the mix.vars required is formed, where for categorical/factor variables a dummy value representation is used.
#'
#' @return A list with
#' \describe{
#' \item{x.link }{A matrix with \eqn{n^2} rows and length(link.vars) columns, detailing the differences in covariates for all pairs of nodes.}
#' \item{x.mix }{A matrix with n rows and number of columns equal to the number of variables detailed in mix.vars, where dummy variable representations will be used for categorical.factor covariates.}
#' }
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
#' @examples data(us.twitter.covariates)
#' link.vars = c(1)
#' mix.vars = c(1,5,7,8)
#' res = formatting.covars(us.twitter.covariates, link.vars, mix.vars, nrow(us.twitter.covariates))
#' dim(res$x.link)
#' dim(res$x.mix)
#' @export
formatting.covars <-
function(covars, link.vars, mix.vars, n)
{
  # covariate data for link regression. First column is 1.
  covars.link = as.data.frame(covars[,link.vars], nrow=n)
  x.link = as.data.frame(matrix(1, n^2, length(link.vars)))
  colnames(x.link) = colnames(covars.link)
  if(length(link.vars) > 1){
    for(j in 2:ncol(covars.link))
    {
      if(is.factor(covars.link[,j]))
      {
        x.link[,j] = as.numeric(matrix(as.numeric(covars.link[,j]),n,n)-t(matrix(as.numeric(covars.link[,j]),n,n)) == 0)   # Gives a 1 if the covariates are the same
      }else{
        x.link[,j] = c(matrix(covars.link[,j],n,n)-t(matrix(covars.link[,j],n,n)))
      }
    }}

  # Covariate data for mixing proportions. First column is 1.
  covars.mix = as.data.frame(covars[,mix.vars], nrow=n)
  nl = sapply(covars.mix, nlevels)-1
  nl.total = sum(nl[nl!=-1]) + sum(nl==-1)
  x.mix = as.data.frame(matrix(1,n,nl.total))
  colnames(x.mix)[1] = colnames(covars.mix)[1]
  counter = 2
  if(length(mix.vars) > 1){
    for(j in 2:ncol(covars.mix))
    {
      if(is.factor(covars.mix[,j]))
      {
        for(i in 1:nl[j])
        {
          x.mix[,counter] = as.numeric(covars.mix[,j])==i
          colnames(x.mix)[counter] = colnames(covars.mix)[j]
          counter = counter+1
        }
      }else{
        x.mix[,counter] = covars.mix[,j]
        colnames(x.mix)[counter] = colnames(covars.mix)[j]
        counter = counter+1
      }
    }}

  list(x.link = as.matrix(x.link), x.mix = as.matrix(x.mix))
}
