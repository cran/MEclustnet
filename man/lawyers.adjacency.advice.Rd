\name{lawyers.adjacency.advice}
\alias{lawyers.adjacency.advice}
\docType{data}
\title{
Adjacency matrix detailing the presence or absence of advice links between the `Lazega Lawyers'.
}
\description{
Data on whether or not 71 lawyers in a northeastern American law firm asked each other for advice.
}
\usage{data("lawyers.adjacency.advice")}
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
data(lawyers.adjacency.advice)
str(lawyers.adjacency.advice)
}
\keyword{datasets}
