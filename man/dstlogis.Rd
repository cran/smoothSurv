\name{standardized logistic}
\alias{dstlogis}
\title{
    Density of Standardized Logistic Distribution.
}
\description{
    Density function of the logistic distribution with zero mean and unit variance.
}
\usage{
    dstlogis(x)
}
\arguments{
   \item{x}{
      Vector of quantiles.
   }
}
\details{
   \code{dstlogis(x) = dlogis(x, 0, sqrt(3)/pi)}
}
\value{
   The value of the density.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\seealso{
   \code{\link{dlogis}} for the logistic distribution.
}
\examples{
dstlogis(0)
dstlogis(seq(-3, 3, 0.2))
}
\keyword{distribution}
