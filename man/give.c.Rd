\name{give.c}
\alias{give.c}
\title{
   Work Function for 'smoothSurvReg'
}
\description{
       Give a vector of all 'c' G-spline coefficients 
       satisfying the three constrains (see below)
       if only (g-3) c's are given.
}
\usage{
give.c(knots, sdspline, last.three, c.rest)
}
\arguments{
   \item{knots}{
      A vector of G-spline knots \eqn{\mu}{mu}.
   }
   \item{sdspline}{
      Standard deviation \eqn{\sigma_0}{sigma0} of the basis G-spline.
   }
   \item{last.three}{
      Indeces of the three 'c' G-spline coefficients which are 
      not included in \code{c.rest}. This must be a vector
      of length 3 with three different numbers from
      \code{1:length(knots)}.
   }
   \item{c.rest}{
      A vector of \code{length(knots) - 3} 'c' spline coefficients
      corresponding to knots \code{knots[-last.three]}. 
      They must all lie between zero and one and their sum must be lower 
      or equal to one.
   }
}
\details{
   The functions computes remaining three 'c' G-spline coefficients
   (if it is possible) such that the resulting G-spline satisfies
   the following:
   \deqn{\sum_{j=1}^g c_j = 1,}{c[1] + ... c[g] = 1,}
   \deqn{\sum_{j=1}^g c_j \mu_j = 0,}{c[1]mu[1] + ... c[g]mu[g] = 0,}
   \deqn{\sum_{j=1}^g c_j \mu_j^2 = 1 - \sigma_0^2.}{c[1]mu[1]^2 + ... c[g]mu[g]^2 = 1 - sigma0^2.}
}
\value{
   A vector of all \code{length(knots)} 'c' G-spline coefficients.
}
\note{
WARNING: This function is primarily used inside \code{\link{smoothSurvReg}}.
Consequently, it has very few error checks on its input arguments.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\examples{
knots <- seq(-4, 4, 0.5)
sd0 <- 0.3
ccoef <- find.c(knots, sd0, dist = "dstlogis")

last.three <- c(3, 7, 10)
c.rest <- ccoef[-last.three]
ccoef2 <- give.c(knots, sd0, last.three, c.rest)

print(ccoef)
print(ccoef2)    ## Almost no change
}
\keyword{internal}
\keyword{utilities}
