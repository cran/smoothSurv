\name{extreme value}
\alias{dextreme}
\alias{dstextreme}
\title{
    Density of the Extreme Value Distribution of a Minimum.
}
\description{
    Density function of the extreme value distribution of a minimum
    with location \eqn{\alpha}{alpha} and scale \eqn{\beta}{beta}
    and the density of the standardized version (with zero mean and unit variance).
}
\usage{
dextreme(x, alpha=0, beta=1)
dstextreme(x)
}
\arguments{
   \item{x}{
      Vector of quantiles.
   }
   \item{alpha}{
      Vector of location parameters.
   }
   \item{beta}{
      Vector of scale parameters.
   }
}
\details{
   Extreme value distribution of a minimum with the location \eqn{\alpha}{alpha}
   and the scale \eqn{\beta}{beta} has a density
   \deqn{f(x) = \frac{1}{\beta}\exp\left[\frac{x-\alpha}{\beta}-\exp\left(\frac{x-\alpha}{\beta}\right)\right]}{%
         f(x) = (1/beta)*exp((x-alpha)/beta - exp((x-alpha)/beta))}
   the mean equal to \eqn{\alpha - \beta\;e}{alpha - beta*e}, where \eqn{e}{e} is approximately
   \eqn{0.5772}{0.5772} and the variance equal to \eqn{\beta^2\frac{\pi}{6}}{beta^2 pi^2/6}.
   Its standardized version is obtained with \eqn{\alpha = \frac{\sqrt{6}}{\pi}\;e}{alpha = (sqrt(6)/pi)*e} 
   and \eqn{\beta = \frac{\sqrt{6}}{\pi}}{beta = (sqrt(6)/pi).}
}
\value{
   The value of the density.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.ac.be}
}
\examples{
dextreme(1, (sqrt(6)/pi)*0.5772, sqrt(6)/pi)
dstextreme(1)        ## approximately same result as on the previous row
}
\keyword{distribution}
