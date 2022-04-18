\name{confint.smoothSurvReg}
\alias{confint.smoothSurvReg}
\title{
   Confidence Intervals for Regression Parameters of 'smoothSurvReg' Model
}
\description{
  Computes confidence intervals for one or more regression related
  parameters (regression coefficients, scale parameter or regression
  coefficients in a model for scale) for a 'smoothSurvReg' model.
}
\usage{
\method{confint}{smoothSurvReg}(object, parm, level = 0.95,
    method = c("pseudo-variance", "asymptotic"), \dots)
}
\arguments{
   \item{object}{
      Object of class smoothSurvReg.
   }
   \item{parm}{A specification of which parameters are to be given
     confidence intervals, either a vector of numbers or a vector of
     names. If missing, all parameters are considered.
   }
   \item{level}{The confidence level required.
   }  
   \item{method}{Type of confidence intervals to b calculated. Option
     \dQuote{pseudo-variance} provides confidence intervals derived from
     inverted minus second derivatives of the penalized log-likelihood
     (pseudo-variance matrix), option \dQuote{asymptotic} provides
     confidence intervals derived from the asyptotic covariance matrix
     of the parameter estimates.
   }  
   \item{\dots}{
      Argument included in the function parameters for the compatibility 
      with the generic function.
   }
}
\value{
   A matrix (or vector) with columns giving lower and upper confidence
   limits for each parameter. These will be labelled as (1 - level)/2
   and 1 - (1 - level)/2 in \% (by default 2.5 \% and 97.5 \%).
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\seealso{
   \code{\link{smoothSurvReg}}
}
\keyword{methods}
