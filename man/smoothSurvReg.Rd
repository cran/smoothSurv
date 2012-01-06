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
smoothSurvReg(formula = formula(data), logscale = ~1, 
   data = parent.frame(), subset, na.action = na.fail,
   init.beta, init.logscale, init.c, init.dist = "best",
   update.init = TRUE, aic = TRUE, lambda = exp(2:(-9)),
   model = FALSE, control = smoothSurvReg.control(), \dots)
}
\arguments{
   \item{formula}{
      A formula expression as for other regression models.
      See the documentation for \code{\link{lm}} and \code{\link{formula}} for details.
      Use \code{\link[survival]{Surv}} on the left hand side of the formula.
   }
   \item{logscale}{
      A formula expression to determine a possible dependence of the
      log-scale on covariates.
   }     
   \item{data}{
      Optional data frame in which to interpret the variables occurring in the
      formula.
   }
   \item{subset}{
      Subset of the observations to be used in the fit.
   }
   \item{na.action}{
      Function to be used to handle any NAs in the data. It's default
      value is \code{na.fail}. It is not recommended to change it in the
      case when \code{logscale} depends on covariates.
   }
   \item{init.beta}{
      Optional vector of the initial values of the regression parameter \eqn{\beta}{beta}
      (intercept and regression itself).
   }
   \item{init.logscale}{
      Optional value of the initial value of the parameters that
      determines the log-scale parameter \eqn{\log(\sigma)}{log(sigma)}.
   }
   \item{init.c}{
      Optional vector of the initial values for the G-spline
      coefficients c, all values must lie between 0 and 1
      and must sum up to 1.
   }
   \item{init.dist}{
      A character string specifying the distribution used by \code{\link[survival]{survreg}}
      to find the initial values for parameters (if not given by the user).
      It is assumed to name \code{"best"} or an
      element from \code{\link[survival]{survreg.distributions}}. These include
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
      \code{\link[survival]{survreg}} fit as determined by the parameter
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
   Read the papers referred below.

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
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\examples{
##### EXAMPLE 1:  Common scale
##### ========================
### We generate interval censored data and fit a model with few artificial covariates
set.seed(221913282)
x1 <- rbinom(50, 1, 0.4)                                         ## binary covariate
x2 <- rnorm(50, 180, 10)                                         ## continuous covariate
y1 <- 0.5*x1 - 0.01*x2 + 0.005 *x1*x2 + 1.5*rnorm(50, 0, 1)      ## generate log(T), left limit
t1 <- exp(y1)                                                    ## left limit of the survival time
t2 <- t1 + rgamma(50, 1, 1)                                      ## right limit of the survival time
surv <- Surv(t1, t2, type = "interval2")                         ## survival object

## Fit the model with an interaction
fit1 <- smoothSurvReg(surv ~ x1 * x2, logscale = ~1, info = FALSE, lambda = exp(2:(-1)))

## Print the summary information
summary(fit1, spline = TRUE)

## Plot the fitted error distribution
plot(fit1)

## Plot the fitted error distribution with its components
plot(fit1, components = TRUE)

## Plot the cumulative distribution function corresponding to the error density
survfit(fit1, cdf = TRUE)

## Plot survivor curves for persons with (x1, x2) = (0, 180) and (1, 180)
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
survfit(fit1, cov = cov)

## Plot hazard curves for persons with (x1, x2) = (0, 180) and (1, 180)
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
hazard(fit1, cov = cov)

## Plot densities for persons with (x1, x2) = (0, 180) and (1, 180)
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
fdensity(fit1, cov = cov)

## Compute estimates expectations of survival times for persons with
## (x1, x2) = (0, 180), (1, 180), (0, 190), (1, 190), (0, 200), (1, 200)
## and estimates of a difference of these expectations:
## T(0, 180) - T(1, 180), T(0, 190) - T(1, 190), T(0, 200) - T(1, 200),
cov1 <- matrix(c(0, 180, 0,   0, 190, 0,   0, 200, 0), ncol = 3, byrow = TRUE)
cov2 <- matrix(c(1, 180, 180,   1, 190, 190,   1, 200, 200), ncol = 3, byrow = TRUE)
print(estimTdiff(fit1, cov1 = cov1, cov2 = cov2))


##### EXAMPLE 2:  Scale depends on covariates
##### =======================================
### We generate interval censored data and fit a model with few artificial covariates
set.seed(221913282)
x1 <- rbinom(50, 1, 0.4)                                         ## binary covariate
x2 <- rnorm(50, 180, 10)                                         ## continuous covariate
x3 <- runif(50, 0, 1)                                            ## covariate for the scale parameter
logscale <- 1 + x3
scale <- exp(logscale)
y1 <- 0.5*x1 - 0.01*x2 + 0.005 *x1*x2 + scale*rnorm(50, 0, 1)    ## generate log(T), left limit
t1 <- exp(y1)                                                    ## left limit of the survival time
t2 <- t1 + rgamma(50, 1, 1)                                      ## right limit of the survival time
surv <- Surv(t1, t2, type = "interval2")                         ## survival object

## Fit the model with an interaction
fit2 <- smoothSurvReg(surv ~ x1 * x2, logscale = ~x3, info = FALSE, lambda = exp(2:(-1)))

## Print the summary information
summary(fit2, spline = TRUE)

## Plot the fitted error distribution
plot(fit2)

## Plot the fitted error distribution with its components
plot(fit2, components = TRUE)

## Plot survivor curves for persons with (x1, x2) = (0, 180) and (1, 180)
## x3 = 0.8 and 0.9
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
logscale.cov <- c(0.8, 0.9)
survfit(fit2, cov = cov, logscale.cov = logscale.cov)

## Plot hazard curves for persons with (x1, x2) = (0, 180) and (1, 180)
## x3 = 0.8 and 0.9
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
logscale.cov <- c(0.8, 0.9)
hazard(fit2, cov = cov, logscale.cov=c(0.8, 0.9))

## Plot densities for persons with (x1, x2) = (0, 180) and (1, 180)
## x3 = 0.8 and 0.9
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
logscale.cov <- c(0.8, 0.9)
fdensity(fit2, cov = cov, logscale.cov = logscale.cov)


## More involved examples can be found in script files
## used to perform analyses  and draw pictures 
## presented in above mentioned references.
## These scripts and some additional files can be found as *.tar.gz files
## in the /inst/doc directory of this package.
##
}
\keyword{survival}
\keyword{smooth}
