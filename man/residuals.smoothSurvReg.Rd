\name{residuals.smoothSurvReg}
\alias{residuals.smoothSurvReg}
\title{
   Residuals for Objects of Class 'smoothSurvReg'
}
\description{
   Compute residuals for the fitted model.
}
\usage{
\method{residuals}{smoothSurvReg}(object, \dots)
}
\arguments{
   \item{object}{
      Object of class smoothSurvReg.
   }
   \item{\dots}{
      Argument included in the function parameters for the compatibility 
      with the generic function.
   }
}
\value{
   A dataframe with columns named \code{res}, \code{res2} and \code{censor}
   where the column \code{res2} is included only if there are any interval censored
   observations. Column \code{res} contains all residuals, column \code{res2}
   is applicable only for interval censored observations. Column \code{censor}
   gives the type of censoring (0 for right censoring, 1 for exact observations,
   2 for left censoring and 3 for interval censoring).
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\seealso{
   \code{\link{smoothSurvReg}}
}
\keyword{methods}
