\name{print.smoothSurvReg}
\alias{print.smoothSurvReg}
\alias{summary.smoothSurvReg}
\title{
   Summary and Print for Objects of Class 'smoothSurvReg'
}
\description{
  Print a summary information of the fitted model.

  For \bold{regression coefficients} the following information is given:
  \tabular{lll}{
  \samp{Value} \tab - \tab estimate of the coefficient \cr
  \samp{Std.Error} \tab - \tab estimated standard error based on the
    pseudo-variance estimate (3.1) \cr
    \tab \tab in
         Komárek, Lesaffre and Hilton (2005)
    \cr
  \samp{Std.Error2} \tab - \tab estimated standard error based on the
     asymptotic variance estimate (3.2) \cr
     \tab \tab in
          Komárek, Lesaffre and Hilton (2005)
     \cr
  \samp{Z} \tab - \tab the Wald statistic obtained as \samp{Value} divided by
    \samp{Std.Error} \cr
  \samp{Z2} \tab - \tab the Wald statistic obtained as \samp{Value} divided
    by \samp{Std.Error2} \cr
  \samp{p} \tab - \tab the two-sided P-value based on normality of the statistic
    \samp{Z} \cr
  \samp{p2} \tab - \tab the two-sided P-value based on normality of the statistic
    \samp{Z2} \cr
  }

  Further, we print:
  \tabular{lll}{  
    \samp{Lambda} \tab - \tab the optimal value of the smoothing
    hyperparameter \cr
    \tab \tab divided by the sample size,
    i.e., \eqn{\lambda/n}{lambda/n} in the notation \cr
    \tab \tab of
    Komárek, Lesaffre and Hilton (2005)
    \cr
    \samp{Log(Lambda)} \tab - \tab logarithm of the above \cr
    \samp{df} \tab - \tab effective degrees of freedom of the model, see
    Section 2.2.3 \cr
        \tab \tab of Komárek, Lesaffre and Hilton (2005)
    \cr
    \samp{AIC} \tab - \tab Akaike's information criterion of the model, see
    Section 2.2.3 \cr
      \tab \tab of Komárek, Lesaffre and Hilton (2005)
    \cr
  }  

  With argument \bold{spline} set to \code{TRUE}, analogous table like
  that for the regression coefficients is printed also for the weights of the
  penalized Gaussian mixture (G-spline).           
}
\usage{
\method{print}{smoothSurvReg}(x, spline, digits = min(options()$digits, 4), \dots)
\method{summary}{smoothSurvReg}(object, spline, digits = min(options()$digits, 4), \dots)
}
\arguments{
   \item{x}{Object of class smoothSurvReg.}
   \item{object}{Object of class smoothSurvReg.}
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
\references{
  Komárek, A., Lesaffre, E., and Hilton, J. F. (2005).
  Accelerated failure time model for arbitrarily censored data with smoothed error distribution. 
  \emph{Journal of Computational and Graphical Statistics}, \bold{14},
  726--745.

  Lesaffre, E., Komárek, A., and Declerck, D. (2005).
  An overview of methods for interval-censored data with an emphasis on applications in dentistry.
  \emph{Statistical Methods in Medical Research}, \bold{14}, 
  539--552.  
} 
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\seealso{
   \code{\link{smoothSurvReg}}, \code{\link{print}}, \code{\link{summary}}
}
\keyword{methods}
