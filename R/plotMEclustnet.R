#' Plot latent position network.
#'
#' @description Function to plot the resulting fitted network, using first two dimensions only.
#' @param fit An object storing the output of the function \code{\link{MEclustnet}}.
#' @param Y The n x n binary adjacency matrix, with 0 down the diagonal, that was passed to \code{\link{MEclustnet}}.
#' @param link.vars A vector of the column numbers of the data frame \code{covars} to be included in link probability model. If none are to be included, this argument should be 1.
#' @param mix.vars A vector of the column numbers of the data frame \code{covars} to be included in mixing proportions model. If none are to be included, argument should be 1.
#'
#' @details This function will plot the posterior mean latent location for each node in the network. The colour of each node reflects the posterior modal cluster membership, and the ellipses are 50\% posterior sets illustrating the uncertainty in the latent locations. The grey lines illustrate the observed links between the nodes.
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
#' @importFrom graphics plot segments points
#' @importFrom stats var sd
#' @importFrom ellipse ellipse
#' @examples #################################################################
#' # An example analysing a 2016 Twitter network of US politicians.
#' #################################################################
#' # Number of iterations etc. are set to low values for illustrative purposes.
#' # Longer run times are likely to be required to achieve sufficient mixing.
#'
#' library(latentnet)
#' data(us.twitter.adjacency)
#' data(us.twitter.covariates)
#'
#' link.vars = c(1)
#' mix.vars = c(1,5,7,8)
#'
#' \donttest{fit = MEclustnet(us.twitter.adjacency, us.twitter.covariates,
#'            link.vars, mix.vars, G=4, d=2, itermax = 500, burnin = 50, uphill = 1, thin=10)
#'
#' # Plot the trace plot of the mean of dimension 1 for each cluster.
#' matplot(t(fit$mustore[,1,]), type="l", xlab="Iteration", ylab="Parameter")
#'
#' # Compute posterior summaries
#' summ = summaryMEclustnet(fit, us.twitter.adjacency)
#'
#' plot(summ$zmean, col=summ$Kmode, xlab="Dimension 1", ylab="Dimension 2", pch=summ$Kmode,
#'      main = "Posterior mean latent location for each node.")
#'
#' # Plot the resulting latent space, with uncertainties
#' plotMEclustnet(fit, us.twitter.adjacency, link.vars, mix.vars)
#'
#' # Examine which politicians are in which clusters...
#' clusters = list()
#' for(g in 1:fit$G)
#' {
#'   clusters[[g]] = us.twitter.covariates[summ$Kmode==g,c("name", "party")]
#' }
#' clusters
#'
#' }
#' @export
plotMEclustnet <-
function(fit, Y, link.vars, mix.vars)
{
  n = nrow(Y)

  zmean = apply(fit$zstore,c(1:2),mean,na.rm=TRUE)
  Kmode = as.numeric(apply(fit$Kstore,2,function(v){names(sort(-table(v)))[1]}))

  xlims = c(min(fit$zstore[,1,]), max(fit$zstore[,1,]))
  ylims = c(min(fit$zstore[,2,]), max(fit$zstore[,2,]))
  plot(zmean,type="n",xlim=xlims, ylim=ylims,col=Kmode,xlab="Dimension 1",
       ylab="Dimension 2",
  main = paste("Link covariates:", paste(link.vars, collapse =""), "\n", "Mixing proportion covariates:", paste(mix.vars, collapse ="")))

  mat = matrix(1:n,n,n)
  edges = cbind(mat[Y==1],t(mat)[Y==1])
  for(m in 1:nrow(edges))
  {
    segments(zmean[edges[m,1],1],zmean[edges[m,1],2],zmean[edges[m,2],1],zmean[edges[m,2],2],col="lightgray",lwd=0.75)
  }

  zvar = array(NA,c(2,2,n))
  for(i in 1:n)
  {
    zvar[,,i] = var(t(fit$zstore[i,,]),na.rm=TRUE)
    points(ellipse(zvar[,,i],centre=zmean[i,],level=0.5),type="l",col=Kmode[i])
  }
  points(zmean, col=Kmode, pch=Kmode+1, bg=Kmode, cex=1)
}
