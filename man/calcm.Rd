\name{calcm}
\alias{calcm}

\title{
Totals the number of observations in each cluster.
}
\description{
Update the count of the number of observations in each cluster.
}
\usage{
calcm(m, G, K)
}

\arguments{
  \item{m}{
Vector of length G containing the number of nodes in each cluster.
}
  \item{G}{
The number of clusters in the model being fitted.
}
  \item{K}{
Vector of length n detailing the number of the cluster to which each node belongs.
}
}
\details{
Function used within the \code{\link{MEclustnet}} function.
}
\value{
Vector of length G containing the number of nodes in each cluster.
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
