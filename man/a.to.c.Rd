\name{a.to.c}
\alias{a.to.c}
\title{
   Work Function for 'smoothSurvReg'
}
\description{
   Compute the 'c' G-spline coefficients from 'a' G-spline coefficients.
}
\usage{
a.to.c(acoef)
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
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.ac.be}
}
\seealso{\code{\link{c.to.a}}}
\examples{
ccoef <- c(0.1, 0.2, 0.15, 0.3, 0.25)

### Compute 'a' counterparts
acoef <- c.to.a(ccoef, 1)
print(acoef)

### And back 'c', ccoef2 should be same as ccoef
ccoef2 <- a.to.c(acoef)
print(ccoef2)
}
\keyword{utilities}
