\name{fdensity.smoothSurvReg}
\alias{fdensity.smoothSurvReg}
\alias{fdensity}
\title{
   Density for Objects of Class 'smoothSurvReg'
}
\description{
  Compute and plot density function for given combinations
  of covariates based on the fitted model.
}
\usage{
\method{fdensity}{smoothSurvReg}(x, cov, logscale.cov, time0 = 0, plot = TRUE,
    by, xlim, ylim, xlab = "t", ylab = "f(t)", 
    type = "l", lty, main, sub, legend, bty = "n", cex.legend = 1, \dots)
}
\arguments{
   \item{x}{
      Object of class smoothSurvReg.
   }
   \item{cov}{
      Vector or matrix with covariates values for which the survivor curve/cdf
      is to be computed and plotted. It must be a matrix with as many columns as
      is the number of covariates (interactions included) or the vector of length
      equal to the number of covariates (interactions included). Intercept is not
      to be included in \code{cov}. If \code{cov} is missing a survivor curve
      for the value of a covariate vector equal to zero is plotted. If there is
      only intercept in the model the survivor curve based on the fitted error
      distribution is always plotted.
   }
   \item{logscale.cov}{
      Vector or matrix with covariate values for the expression of log-scale
      (if this depended on covariates).
      It can be omitted in the case that log-scale was common for all observations. 
   }     
   \item{time0}{
      Starting time of the follow-up as used in the model. I.e. the
      model is assumed to be \eqn{\log(T-time0) = x'\beta + 
	\sigma\varepsilon}{log(T-time0) = x'beta + sigma*epsilon} 
   } 
   \item{plot}{
      If \code{TRUE} the plot is directly produced, otherwise only a \code{data.frame}
      with information used for later plotting is returned.
   }
   \item{by}{
      Step for a ploting grid. If \code{missing} it is automatically computed.
   }
   \item{xlim, ylim}{
      Arguments passed to the \code{\link{plot}} function.
   }
   \item{xlab, ylab}{
      Arguments passed to the \code{\link{plot}} function.
   }
   \item{type, lty}{
      Arguments passed to the \code{\link{plot}} function.
   }
   \item{main, sub}{
      Arguments passed to the \code{\link{plot}} function.
   }
   \item{legend, bty}{
      Argument passed to the \code{\link{plot}} function.
    }
   \item{cex.legend}{argument passed to \code{cex} argument of the
     \code{\link{legend}} function.}
   \item{\dots}{
      Arguments passed to the \code{\link{plot}} function.
   }
}
\value{
   A dataframe with columns named \code{x} and \code{y} where \code{x} gives the grid
   and \code{y} the values of the density function at that grid.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\seealso{
   \code{\link{smoothSurvReg}}, \code{\link{plot}}
}
\keyword{methods}
