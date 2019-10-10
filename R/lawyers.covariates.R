#' A matrix of covariates of the `Lazega Lawyers'.
#'
#' Covariates on each of 71 lawyers in a northeastern American law firm.
#' Note the first column is a column of 1's.
#'
#' @format A data frame with 71 observations on the following 8 variables.
#' \describe{
#' \item{\code{Intercept}}{a column of 1s should always be the first column.}
#' \item{\code{Seniority}}{a factor with levels \code{1} = partner, \code{2} = associate.}
#' \item{\code{Gender}}{a factor with \code{1} = male, \code{2} = female.}
#' \item{\code{Office}}{a factor with levels \code{1} = Boston, \code{2} = Hartford and \code{3} = Providence}
#' \item{\code{Years}}{a numeric vector detailing years with the firm.}
#' \item{\code{Age}}{a numeric vector detailing the age of each lawyer.}
#' \item{\code{Practice}}{a factor with levels \code{1} = litigation and \code{2} = corporate.}
#' \item{\code{School}}{a factor with levels \code{1} = Harvard or Yale, \code{2} = University of Connecticut and \code{3} = Other.}
#' }
#' @source E. Lazega, The Collegial Phenomenon: The Social Mechanisms of Cooperation
#'  Among Peers in a Corporate Law Partnership, Oxford University Press, Oxford, England, 2001.
"lawyers.covariates"
