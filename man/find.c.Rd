\name{find.c}
\alias{find.c}
\title{
   Work Function for 'smoothSurvReg'
}
\description{
       Find mixture proportions that approximate
       given distribution by a G-spline mixture.
}
\usage{
find.c(knots, sdspline, dist = "dnorm")
}
\arguments{
   \item{knots}{
      A vector of G-spline knots \eqn{\mu}{mu}.
   }
   \item{sdspline}{
      Standard deviation \eqn{\sigma_0}{sigma0} of the basis G-spline.
   }
   \item{dist}{
      A character string specifying the function used to compute
      a density of the distribution
      you want to approximate.
   }
}
\details{
   The function finds the G-spline coefficients that approximates
   a density given by \code{dist} in such sense that the value
   of the G-spline is exactly equal to the value of that density
   in \code{knots}.
}
\value{
   Either the vector of G-spline 'c' coefficients or \code{NULL}
   if there are problems to find them.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\examples{
knots <- seq(-4, 4, 0.5)
sd0 <- 0.3
ccoef <- find.c(knots, sd0, dist = "dstlogis")

### We plot the approximation together with the truth
###
grid <- seq(-4, 4, 0.05)
truth <- dstlogis(grid)

### Following lines compute the values of the approximation
grid.big <- matrix(grid, nrow = length(grid), ncol = length(knots))
knots.big <- matrix(knots, nrow = length(grid), ncol = length(knots), byrow = TRUE)
normals <- dnorm(grid.big, mean = knots.big, sd = sd0)
approx <- normals \%*\% ccoef

### Plot it
plot(grid, approx, type = "l", xlab = "y", ylab = "f(y)", bty = "n")
lines(grid, truth, lty = 2)
legend(-4, 0.35, c("approx", "truth"), lty = 1:2, bty = "n")
}
\keyword{internal}
\keyword{utilities}
