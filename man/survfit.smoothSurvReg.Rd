\name{survfit.smoothSurvReg}
\alias{survfit.smoothSurvReg}
\title{
   Survivor Curves for Objects of Class 'smoothSurvReg'
}
\description{
   Compute and plot survivor function/cumulative distribution function
   for the fitted model.
}
\usage{
survfit.smoothSurvReg(formula, cov = NULL, plot = TRUE, cdf = FALSE, by = NULL,
   xlim = NULL, ylim = c(0, 1), xlab = "t", ylab = expression(paste("S(","t",")", sep = "")),
   type = "l", lty = NULL, main = NULL, sub = NULL, legend = NULL, bty = "n", \dots)
}
\arguments{
   \item{formula}{
      Object of class smoothSurvReg.
   }
   \item{cov}{
      Vector or matrix with covariates values for which the survivor curve/cdf
      is to be computed and plotted. It must be a matrix with as many columns as
      is the number of covariates (interactions included) or the vector of length
      equal to the number of covariates (interactions included). Intercept is not
      to be included in \code{cov}. If \code{cov == NULL} a survivor curve
      for the value of a covariate vector equal to zero is plotted. If there is
      only intercept in the model the survivor curve based on the fiiied error
      distribution is always plotted.
   }
   \item{plot}{
      If \code{TRUE} the plot is directly produced, otherwise only a \code{data.frame}
      with information used for later plotting is returned.
   }
   \item{cdf}{
     If \code{TRUE} cumulative distribution function is plotted instead of
     the survivor function.
   }
   \item{by}{
      Step for a ploting grid. If \code{NULL} it is automatically computed.
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
   \item{\dots}{
      Arguments passed to the \code{\link{plot}} function.
   }
}
\value{
   A dataframe with columns named \code{x} and \code{y} where \code{x} gives the grid
   and \code{y} the values of the survivor/cum. distribution function at that grid.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.ac.be}
}
\seealso{
   \code{\link{smoothSurvReg}}, \code{\link{plot}}
}
\keyword{methods}
