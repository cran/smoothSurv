\name{plot.smoothSurvReg}
\alias{plot.smoothSurvReg}
\title{
   Plot Objects of Class 'smoothSurvReg'
}
\description{
   Plot the fitted error distribution.
}
\usage{
\method{plot}{smoothSurvReg}(x, plot = TRUE, resid = TRUE, knots = TRUE,
   compare = TRUE, components = FALSE, standard = TRUE,
   by, toler.c = 1e-5,
   xlim, ylim,
   xlab = expression(epsilon), ylab = expression(paste("f(",epsilon,")", sep = "")),
   type = "l", lty = 1, main, sub, bty = "n", \dots)
}
\arguments{
   \item{x}{
      Object of class smoothSurvReg.
   }
   \item{plot}{
      If \code{TRUE} the plot is directly produced, otherwise only a \code{data.frame}
      with information used for later plotting is returned.
   }
   \item{resid}{
      If \code{resid & plot} residuals are added
      to the plot in the following way. Residuals based on exact observations
      are plotted by \code{pch = 3} (+), residuals based on right censored
      observations using \code{pch = 4} (x), residuals based on left censored
      observations using \code{pch = 2} (traingle) and midpoints of
      residuals based on interval censored observations using
      \code{pch = 5} (diamond).
   }
   \item{knots}{
      If \code{knots & plot} bullets indicating (transformed)
      knots corresponding to 'c' G-spline coefficients higher than \code{toler.c}
      are added to the x-axis of the plot.
   }
   \item{compare}{
      If \code{compare & !components & standard & plot} plots of the three
      parametric distributions (st. normal, st. logistic and st. minimum extreme value)
      are added to the plot for comparison.
   }
   \item{components}{
      If \code{components & standard & plot} dashed lines with weighted
      G-spline components are added to the plot.
   }
   \item{standard}{
      If \code{standard} data for plotting the fitted error distribution
      are created, otherwise data for plotting the distribution
      of \eqn{\alpha + \sigma\varepsilon}{alpha + sigma*epsilon} are created.
      In the case that \eqn{\log(\sigma)}{log(sigma)} depends on covariates
      and \code{standard} is FALSE, all covariates equal to zero are used
      to compute \eqn{\log(\sigma)}{log(sigma)}.
   }
   \item{by}{
      Step for a ploting grid. If \code{NULL} it is automatically computed.
   }
   \item{toler.c}{
      Tolerance to indicate zero 'c' G-spline coefficients used if \code{knots == TRUE}.
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
   \item{bty}{
      Argument passed to the \code{\link{plot}} function.
   }
   \item{\dots}{
      Arguments passed to the \code{\link{plot}} function.
   }
}
\value{
   A dataframe with columns named \code{x} and \code{y} where \code{x} gives the grid
   and \code{y} the values of the density at that grid.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.ac.be}
}
\seealso{
   \code{\link{smoothSurvReg}}, \code{\link{plot}}
}
\keyword{methods}
