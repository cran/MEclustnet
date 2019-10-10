#' Account for invariance of configurations.
#'
#'This function accounts for the fact that configurations in the latent space are invariant to rotations, reflections and translations.
#'
#' @param z An n x d matrix of latent locations in the d dimensional space for each of n nodes.
#' @param zMAP The \emph{maximum a posteriori} configuration of latent locations used as the template to which all sampled configurations are mapped.
#'
#' @return The transformed version of the input configuration z that best matches zMAP.
#' @details Procrustean rotations, reflections and translations (note: NOT dilations) are employed to best match z to zMAP.
#' @seealso \code{\link{MEclustnet}}
#' @references Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
#' @importFrom vegan procrustes
invariant <-
function(z, zMAP)
{
  procrustes(zMAP, z, scale = FALSE)$Yrot
}
