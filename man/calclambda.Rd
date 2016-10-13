\name{calclambda}
\alias{calclambda}

\title{
Compute mixing proportions
}
\description{
Function to compute the each observation's mixing proportions which are modeled as a logistic function of their covariates.
}
\usage{
calclambda(tau, x.mix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau}{
A matrix of logistic regression coefficients, with G rows and number of columns equal to the number of covariates in the mixing proportions model plus 1, for the intercept.
}
  \item{x.mix}{
A matrix of covariates in the mixing proportions model (including dummy variables for any factor covariates), with a column of 1's appended at the front.
}
}
\details{
Function used within the \code{\link{MEclustnet}} function.
}
\value{
An n x G matrix of mixing proportions.
}
\references{
Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
}
\author{
Isobel Claire Gormley <claire.gormley@ucd.ie>
}

\seealso{
\code{\link{MEclustnet}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
