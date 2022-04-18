\name{print.estimTdiff}
\alias{print.estimTdiff}
\title{
   Print for Objects of Class 'estimTdiff'
}
\description{
   Print a summary information of the estimates and tests for expected
   values of survival times based on a regression.
}
\usage{
\method{print}{estimTdiff}(x, digits = min(options()$digits, 4), \dots)
}
\arguments{
   \item{x}{
      Object of class estimTdiff.
   }
   \item{digits}{
      Controls the number of digits to print when printing numeric
      values.  It is a suggestion only.  Valid values are 1...22.
   }
   \item{\dots}{
      Further arguments passed to or from other methods.
   }
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\seealso{
   \code{\link{smoothSurvReg}}, \code{\link{print}}
}
\keyword{methods}
