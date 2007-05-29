\name{estimTdiff}
\alias{estimTdiff}
\alias{estimTdiff.smoothSurvReg}
\title{
  Estimate expectation of survival times and their difference from the
  results given by survival regression function
}
\description{
  Estimate expectation of survival times and their difference from the
  results given by survival regression function
}
\usage{
estimTdiff(x, \dots)
\method{estimTdiff}{smoothSurvReg}(x, cov1, cov2, logscale.cov1, logscale.cov2, time0 = 0, \dots)
}
\arguments{
   \item{x}{
      Object of an appropriate class.
   }
   \item{cov1}{
      Vector or matrix with covariates values for which the expectations
      of the first survival time are to be computed. It must be a matrix with as many 
      columns as is the number of covariates (interactions included, intercept excluded) 
      or the vector of length
      equal to the number of covariates (interactions included,
      intercept excluded). If matrix is supplied then is assumed that
      each row of this matrix gives one covariate combination for the
      first survival time.
      Intercept is not to be included in \code{cov1}. If \code{cov1} is missing an
      expectation of a survivor time
      for the value of a covariate vector equal to zero is computed. If there is
      only intercept in the model, this parameter must be always missing.
    }
    \item{cov2}{
      Vector or matrix with covariate values for which the expectations
      of the second survival time are to be computed. It must be of same
      size as \code{cov1}.
    }          
    \item{logscale.cov1}{
      Vector or matrix with covariate values for the expression of log-scale
      (if this depended on covariates) for the first survival time. 
      It can be omitted in the case that log-scale was common for all observations. 
    } 
    \item{logscale.cov2}{
      Vector or matrix with covariate values for the expression of log-scale
      (if this depended on covariates) for the second survival time. 
      It can be omitted in the case that log-scale was common for all observations. 
    }
    \item{time0}{
      Starting time of the follow-up as used in the model. I.e. the
      model is assumed to be \eqn{\log(T-time0) = x'\beta + 
	\sigma\varepsilon}{log(T-time0) = x'beta + sigma*epsilon} 
    }   
    \item{\dots}{who knows
    }
}
\value{
   A \code{data.frame} with columns named ``ET1'', ``sd.ET1'', ``ET2'',
   ``sd.ET2'', ``diffT'', ``sd.diffT'' giving the estimates of expected
   values of the survival times for covariate values given in rows of
   \code{cov1} and \code{logscale.cov1}, 
   their standard errors, estimates of expected values of
   survival times for covariate values given in rows of \code{cov2} 
   and \code{logscale.cov2},
   their standard errors and estimates of a difference of expected
   values of survival times for covariate values given in rows of
   \code{cov1}, \code{logscale.cov1} and \code{cov2}, \code{logscale.cov2} 
   and their standard errors.
}   
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek@med.kuleuven.be}
}
\seealso{
   \code{\link{smoothSurvReg}}
}
\keyword{methods}
