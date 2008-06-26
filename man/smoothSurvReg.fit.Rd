\name{smoothSurvReg.fit}
\alias{smoothSurvReg.fit}
\title{
   Work Function to Fit the Model Using 'smoothSurvReg'
}
\description{
   Fit the survival regression model with smoothed error distribution.
   This function is not to be called by the user.
}
\usage{
smoothSurvReg.fit(x, z, y, offset = NULL, correctlik, init, controlvals, common.logscale)
}
\arguments{
   \item{x}{
       A covariate matrix.
   }
   \item{z}{
       A covariate matrix for log(scale).
   }     
   \item{y}{
       A two- or three- column matrix with response. The last column
       indicate the type of censoring. See \code{\link[survival]{Surv}} for details.
   }
   \item{offset}{
       A vector with possible offset term.
   }
   \item{correctlik}{
       Correction term to the likelihood due to the log transformation of the response.
   }
   \item{init}{
       A list with components \code{beta}, \code{scale}, \code{ccoef}
       giving the initial values for the fitting process.
   }
   \item{controlvals}{
       A list returned by \code{\link{smoothSurvReg.control}}.
   }
   \item{common.logscale}{
       Indicator (TRUE/FALSE) indicating whether the log-scale is same for all observations 
       or whether it depends on covariates.
   }
}
\value{
   A list with components \code{regres, spline, loglik, aic,
     degree.smooth, var, var2, dCdC,
   iter, estimated, warning, fail, H, I, G, U}.
}
\note{
WARNING: Most users will call the higher level routine smoothSurvReg.
Consequently, this function has very few error checks on its input arguments.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}
\keyword{survival}
\keyword{smooth}
