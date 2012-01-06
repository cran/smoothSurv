\name{MP.pseudoinv}
\alias{MP.pseudoinv}
\title{
   Work Function for 'smoothSurvReg', currently nowhere used.
}
\description{
   Function to compute a Moore-Penrose pseudoinverse
   of a symmetric matrix using eigen values decomposition.
}
\usage{
MP.pseudoinv(x, toler = 1e-7)
}
\arguments{
   \item{x}{
      A symmetric matrix (it is not checked).
   }
   \item{toler}{
      Under the eigen values decomposition is computed, all
      eigen values smaller in absolute value than \code{toler}
      are assumed to be zero.
   }
}
\value{
   A matrix with the Moore-Penrose pseudoinverse of \code{x}.
}
\note{
WARNING: This function was primarily used inside \code{\link{smoothSurvReg}}.
Consequently, it has very few error checks on its input arguments.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\keyword{internal}
\keyword{utilities}
