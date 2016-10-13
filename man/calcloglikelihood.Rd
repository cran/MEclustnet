\name{calcloglikelihood}
\alias{calcloglikelihood}

\title{
Calculate the log likelihood function of the data.
}
\description{
This function calculates the log likelihood function of the data.
}
\usage{
calcloglikelihood(pis, y)
}

\arguments{
  \item{pis}{
Vector of link probabilities.
}
  \item{y}{
Vector version of the adjacency matrix, with the diagonal removed.
}
}
\details{
Function used within other functions.
}
\value{
The value of the log likelihood function.
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
