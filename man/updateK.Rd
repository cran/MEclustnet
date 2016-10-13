\name{updateK}
\alias{updateK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update the cluster membership vector.
}
\description{
A Gibbs update step for K, the cluster membership vector.
}
\usage{
updateK(G, K, z, mu, sigma2, Id, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
The number of clusters being fitted.
}
  \item{K}{
The cluster membership vector.
}
  \item{z}{
The n x d matrix of latent locations.
}
  \item{mu}{
The G x d matrix of cluster means.
}
  \item{sigma2}{
The G vector of cluster covariances.
}
  \item{Id}{
An identity matrix of dimension d.
}
  \item{lambda}{
The n x G matrix of mixing proportions.
}
}

\value{
The cluster membership vector.
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
