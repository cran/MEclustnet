\name{updatemu}
\alias{updatemu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update the mean of each cluster.
}
\description{
A Gibbs step to update the mean of each cluster.
}
\usage{
updatemu(G, z, K, m, sigma2, omega2, Id, mu, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
The number of clusters being fitted.
}
  \item{z}{
The n x d matrix of latent locations.
}
  \item{K}{
The cluster membership vector.
}
  \item{m}{
Vector of length G containing the number of nodes in each cluster.
}
  \item{sigma2}{
The covariance of each cluster.
}
  \item{omega2}{
Covariance of the multivariate normal prior distribution on the means. Note this is a scalar value, as the prior covariance is diagonal.
}
  \item{Id}{
A d x d identity matrix.
}
  \item{mu}{
The G x d matrix of cluster means.
}
  \item{d}{
The dimension of the latent space.
}
}

\value{
The G x d matrix of cluster means.
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
