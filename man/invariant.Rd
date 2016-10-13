\name{invariant}
\alias{invariant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Account for invariance of configurations.
}
\description{
This function accounts for the fact that configurations in the latent space are invariant to rotations, reflections and translations.
}
\usage{
invariant(z, zMAP)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
An n x d matrix of latent locations in the d dimensional space for each of n nodes.
}
  \item{zMAP}{
The \emph{maximum a posteriori} configuration of latent locations used as the template to which all sampled configurations are mapped.
}
}
\details{
Procrustean rotations, reflections and translations (note: NOT dilations) are employed to best match z to zMAP.
}
\value{
The transformed version of the input configuration z that best matches zMAP.
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

\keyword{ internal }% use one of  RShowDoc("KEYWORDS")

