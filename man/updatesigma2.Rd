\name{updatesigma2}
\alias{updatesigma2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update variances in each cluster.
}
\description{
A Gibbs step to update variances in each cluster.
}
\usage{
updatesigma2(G, alpha, m, d, sigma02, z, K, mu, sigma2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
The number of clusters being fitted.
}
  \item{alpha}{
Degrees of freedom of the scaled inverse Chi squared prior distribution on the cluster variances.
}
  \item{m}{
Vector of length G containing the number of nodes in each cluster.
}
  \item{d}{
Dimension of the latent space.
}
  \item{sigma02}{
Scaled factor of the scaled inverse Chi squared prior distribution on the cluster variances.
}
  \item{z}{
The n x d matrix of latent locations.
}
  \item{K}{
The cluster membership vector.
}
  \item{mu}{
The G x d matrix of cluster means.
}
  \item{sigma2}{
The G vector of cluster variances.
}
}

\value{
The G vector of cluster variances.
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
