\name{us.twitter.adjacency}
\alias{us.twitter.adjacency}
\docType{data}
\title{
Directed adjacency matrix detailing the presence or absence of Twitter friend/follower links between US politicians.
}
\description{
Network data on whether or not 69 US politicians are friends/followers on Twitter.
}
\usage{data("us.twitter.adjacency")}
\format{
  A 69 x 69 binary matrix, with 0 down the diagonal.
}

\source{
With thanks to Dr. Derek Greene. School of Computer Science, University College Dublin.
}

\examples{
data(us.twitter.adjacency)
str(us.twitter.adjacency)
}
\keyword{datasets}
