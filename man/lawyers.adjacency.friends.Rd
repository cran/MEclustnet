\name{lawyers.adjacency.friends}
\alias{lawyers.adjacency.friends}
\docType{data}
\title{
Adjacency matrix detailing the presence or absence of friendship links between the `Lazega Lawyers'.
}
\description{
Data on whether or not 71 lawyers in a northeastern American law firm are friends outside of work.
}
\usage{data("lawyers.adjacency.friends")}
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
data(lawyers.adjacency.friends)
str(lawyers.adjacency.friends) 
}
\keyword{datasets}
