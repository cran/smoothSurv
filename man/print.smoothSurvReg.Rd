\name{print.smoothSurvReg}
\alias{print.smoothSurvReg}
\alias{summary.smoothSurvReg}
\title{
   Summary and Print for Objects of Class 'smoothSurvReg'
}
\description{
   Print a summary information of the fitted model.
}
\usage{
\method{print}{smoothSurvReg}(x, spline, digits = min(options()$digits, 4), \dots)
\method{summary}{smoothSurvReg}(object, spline, digits = min(options()$digits, 4), \dots)
}
\arguments{
   \item{x, object}{
      Object of class smoothSurvReg.
   }
   \item{spline}{
      TRUE/FALSE. If TRUE an information on fitted G-spline is printed.
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
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\seealso{
   \code{\link{smoothSurvReg}}, \code{\link{print}}, \code{\link{summary}}
}
\keyword{methods}
