\name{smoothSurvReg}
\alias{smoothSurvReg}
\title{
   Regression for a Survival Model with Smoothed Error Distribution
}
\description{
   Regression for a survival model. These are all
   time-transformed location models, with the most useful case being the
   accelerated failure models that use a log transformation.
   Error distribution is assumed to be a mixture of G-splines.
   Parameters are estimated by the penalized maximum likelihood
   method.
}
\usage{
smoothSurvReg(formula = formula(data), data = parent.frame(), subset, na.action,
   init.beta = NULL, init.scale = NULL, init.c = NULL, init.dist = "best",
   update.init = TRUE, aic = TRUE, lambda = exp(2:(-9)),
   model = FALSE, control = smoothSurvReg.control(), \dots)
}
\arguments{
   \item{formula}{
      A formula expression as for other regression models.
      See the documentation for \code{\link{lm}} and \code{\link{formula}} for details.
      Use \code{\link{Surv}} on the left hand side of the formula.
   }
   \item{data}{
      Optional data frame in which to interpret the variables occurring in the
      formula.
   }
   \item{subset}{
      Subset of the observations to be used in the fit.
   }
   \item{na.action}{
      Function to be used to handle any NAs in the data.
   }
   \item{init.beta}{
      Optional vector of the initial values of the regression parameter \eqn{\beta}{beta}
      (intercept and regression itself).
   }
   \item{init.scale}{
      Optional value of the initial value of the scale parameter \eqn{\sigma}{sigma}.
   }
   \item{init.c}{
      Optional vector of the initial values for the G-spline
      coefficients c, all values must lie between 0 and 1
      and must sum up to 1.
   }
   \item{init.dist}{
      A character string specifying the distribution used by \code{\link{survreg}}
      to find the initial values for parameters (if not given by the user).
      It is assumed to name \code{"best"} or an
      element from \code{\link{survreg.distributions}}. These include
      \code{"weibull"}, \code{"exponential"}, \code{"gaussian"},
      \code{"logistic"}, \code{"lognormal"} and \code{"loglogistic"}.
      If "best" is specified one of \code{"lognormal"},
      \code{"weibull"} and \code{"loglogistic"} giving the highest
      likelihood is used.
   }
   \item{update.init}{
      If TRUE, the initial values are updated during the grid search for the lambda
      parameter giving the optimal AIC. Otherwise, fits with all lambdas during the
      grid search start with same initials determine at the beginning either
      from the values of \code{init.beta, init.scale, init.c} or from the initial
      \code{\link{survreg}} fit as determined by the parameter
      \code{init.dist}.
   }
   \item{aic}{
      If TRUE the optimal value of the tuning parameter \eqn{\lambda}{lambda}
      is determined via a grid search through the values specified
      by the parameter \code{lambda}. If FALSE, only the model with
      \eqn{\lambda = }{lambda = } \code{lambda[1]} is fitted.
   }
   \item{lambda}{
      A grid of values of the tuning parameter \eqn{\lambda}{lambda} searched
      for the optimal value if \code{aic} = TRUE.
   }
   \item{model}{
      If TRUE, the model frame is returned.
   }
   \item{control}{
      A list of control values, in the format producted by \code{\link{smoothSurvReg.control}}.
   }
   \item{...}{
      Other arguments which will be passed to \code{\link{smoothSurvReg.control}}. See its help
      page for more options to control the fit and for the possibility to fix some values
      and not to estimate them.
   }
}
\details{
   Read the paper paper.pdf found in /inst/doc directory of this library.

   There is a slight difference in the definition of the penalty used by the R function compared
   to what is written in the paper.
   The penalized log-likelihood given in the paper has a form
   \deqn{\ell_P(\theta) = \ell(\theta) - \frac{\lambda}{2}\sum_{j=m+1}^g(\Delta^m a_j)^2,}{%
        l_P(theta) = l(theta) - (lambda/2) * sum[j in (m+1):g] (Delta^m a[j])^2,}
   while the penalized log-likelihood used in the R function multiplies the tuning parameter
   \eqn{\lambda}{lambda} given by \code{lambda} by a sample size \eqn{n}{n} to keep default values
   more or less useful for samples of different sizes. So that the penalized log-likelihood
   which is maximized by the R function has the form
   \deqn{\ell_P(\theta) = \ell(\theta) - \frac{\lambda\cdot n}{2}\sum_{j=m+1}^g(\Delta^m a_j)^2.}{%
        l_P(theta) = l(theta) - ((lambda*n)/2) * sum[j in (m+1):g] (Delta^m a[j])^2.}
}
\value{
   An object of class \code{smoothSurvReg} is returned.
   See \code{\link{smoothSurvReg.object}} for details.
}
\references{
  \eqn{\mbox{Kom\'arek, A., Lesaffre, E., and Hilton, J. F. (2004)}}{Komarek, A., Lesaffre, E., and Hilton, J. F. (2004)}
  Accelerated Failure Time Model for Arbitrarily censored Data
  with Smoothed Error Distribution.
  \emph{?????}, \bold{??}, ???--???.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.ac.be}
}
\examples{
### We generate interval censored data and fit a model with few artificial covariates
x1 <- rbinom(50, 1, 0.4)                                         ## binary covariate
x2 <- rnorm(50, 180, 10)                                         ## continuous covariate
y1 <- 0.5*x1 - 0.01*x2 + 0.005 *x1*x2 + 1.5*rnorm(50, 0, 1)      ## generate log(T), left limit
t1 <- exp(y1)                                                    ## left limit of the survival time
t2 <- t1 + rgamma(50, 1, 1)                                      ## right limit of the survival time
surv <- Surv(t1, t2, type = "interval2")                         ## survival object

## Fit the model with an interaction
fit1 <- smoothSurvReg(surv ~ x1 * x2, info = FALSE, lambda = exp(2:(-1)))

## Print the summary information
summary(fit1, spline = TRUE)

## Plot the fitted error distribution
plot(fit1)

## Plot the fitted error distribution with its components
plot(fit1, components = TRUE)

## Plot the cumulative distribution function corresponding to the error density
survfit(fit1, cdf = TRUE)

## Plot survivor curves for persons with (x1, x2) = (0, 180) and (1, 180)
survfit(fit1, matrix(c(0, 1, 180, 180, 0, 180), ncol = 3))

## More involved examples can be found in script files
## used to perform analyses  and draw pictures 
## presented in /inst/doc/paper.pdf of this library.
## These script files together with data for Signal Tandmobiel
## example (section 6.1 of paper.pdf) can be obtained from
## arnost.komarek@med.kuleuven.ac.be
}
\keyword{survival}
\keyword{smooth}
