\name{us.twitter.covariates}
\alias{us.twitter.covariates}
\docType{data}
\title{
A matrix of covariates of the US politicians.
}
\description{
Covariates on each of 69 US policitians. Note the first column is a column of 1's.
}
\usage{data("us.twitter.covariates")}
\format{
  A data frame with 69 observations on the following 8 variables.
  \describe{
    \item{\samp{1}}{aa column of 1s should always be the first column.}
    \item{\code{twitter_id}}{Twitter number.}
    \item{\code{twitter_name}}{Twitter name.}
    \item{\code{name}}{Actual name.}
    \item{\code{party}}{a factor with levels \code{Democrat} \code{Republican}}
    \item{\code{location}}{a factor with levels detailing location.}
    \item{\code{role}}{a factor with levels \code{Candidate}, \code{Representative} and \code{Senator}}
    \item{\code{gender}}{a factor with levels \code{Female} and \code{Male}}
  }
}
\source{
With thanks to Dr. Derek Greene. School of Computer Science, University College Dublin.
}
\references{
Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
}
\examples{
data(us.twitter.covariates)
str(us.twitter.covariates)
}
\keyword{datasets}
