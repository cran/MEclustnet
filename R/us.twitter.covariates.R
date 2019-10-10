#' A matrix of covariates of the US politicians.
#'
#' Covariates on each of 69 US politicians. Note the
#' first column is a column of 1's.
#'
#' @format A data frame with 69 observations on the following 8 variables.
#' \describe{
#'   \item{\samp{1}}{a column of 1s should always be the first column.}
#' \item{\code{twitter_id}}{Twitter number.}
#' \item{\code{twitter_name}}{Twitter name.}
#' \item{\code{name}}{Actual name.}
#' \item{\code{party}}{a factor with levels \code{Democrat} \code{Republican}}
#' \item{\code{location}}{a factor with levels detailing location.}
#' \item{\code{role}}{a factor with levels \code{Candidate}, \code{Representative} and \code{Senator}}
#' \item{\code{gender}}{a factor with levels \code{Female} and \code{Male}}
#' }
#' @source E. Lazega, The Collegial Phenomenon: The Social Mechanisms of Cooperation
#'  Among Peers in a Corporate Law Partnership, Oxford University Press, Oxford, England, 2001.
"us.twitter.covariates"
