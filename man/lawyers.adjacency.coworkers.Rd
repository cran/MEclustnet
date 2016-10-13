\name{lawyers.adjacency.coworkers}
\alias{lawyers.adjacency.coworkers}
\docType{data}
\title{
Adjacency matrix detailing the presence or absence of coworker links between the `Lazega Lawyers'.
}
\description{
Data on whether or not 71 lawyers in a northeastern American law firm work with each other.
}
\usage{data("lawyers.adjacency.coworkers")}
\format{
  A 71 x 71 binary matrix, with 0 down the diagonal.
}
\source{
E. Lazega, The Collegial Phenomenon: The Social Mechanisms of Cooperation Among Peers in a Corporate Law Partnership, Oxford University Press, Oxford, England, 2001.
}
\references{
Isobel Claire Gormley and Thomas Brendan Murphy. (2010) A Mixture of Experts Latent Position Cluster Model for Social Network Data. Statistical Methodology, 7 (3), pp.385-405.
}
\examples{
data(lawyers.adjacency.coworkers)
str(lawyers.adjacency.coworkers)
}
\keyword{datasets}
