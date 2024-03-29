\name{a2c}
\alias{a2c}
\title{
   Work Function for 'smoothSurvReg'
}
\description{
   Compute the 'c' G-spline coefficients from 'a' G-spline coefficients.
}
\usage{
a2c(acoef)
}
\arguments{
   \item{acoef}{
       A vector of 'a' G-spline coefficients.
   }
}
\details{
   'c' and 'a' G-spline coefficients are related by the expression
   \deqn{c_j = \frac{\exp(a_j)}{\sum_{l=1}^{g}\exp(a_l)}, j = 1,\dots, g,}{%
         c[j] = exp(a[j])/[exp(a[1]) + ... + exp(a[g])], j = 1, ..., g,}
   where \eqn{a_k = 0,}{a[k] = 0,} k = \code{which.zero}.
   This function transforms the 'c' coefficients into their 'a' counterparts.

}
\value{
   A vector of same length as \code{acoef} with 'c' G-spline coefficients.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\seealso{\code{\link{c2a}}}
\examples{
ccoef <- c(0.1, 0.2, 0.15, 0.3, 0.25)

### Compute 'a' counterparts
acoef <- c2a(ccoef, 1)
print(acoef)

### And back 'c', ccoef2 should be same as ccoef
ccoef2 <- a2c(acoef)
print(ccoef2)
}
\keyword{internal}
\keyword{utilities}
