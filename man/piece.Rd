\name{piece}
\alias{piece}
\title{
    Left Continuous Piecewise Constant Function with a Finite Support.
}
\description{
    Function to evaluate a left continuous piecewise constant function with a finite support.
}
\usage{
    piece(x, breaks, values)
}
\arguments{
   \item{x}{
      Vector of values where the piecewise constant function should be evaluated.
   }
   \item{breaks}{
      Vector of sorted breakpoints of the piecewise constant function.
   }
   \item{values}{
      Values of the piecewise constant function. It takes the value
      \code{value[i]} on the interval \code{(breaks[i], breaks[i+1]]}. 
      The function is assumed to be zero outside its range specified
      as \code{(breaks[1], breaks[length(breaks)]]}.
      The length of the vector \code{values} must be equal to
      \code{length(breaks) - 1}
   }
}
\value{
   The value of the piecewise constant function.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\examples{
my.breaks <- c(-2, 1.5, 4, 7)
my.values <- c(0.5, 0.9, -2)
grid <- seq(-3, 8, by = 0.25)
piece(grid, my.breaks, my.values)
}
\keyword{utilities}
