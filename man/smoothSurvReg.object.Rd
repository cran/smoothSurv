\name{smoothSurvReg.object}
\alias{smoothSurvReg.object}
\title{
   Smoothed Survival Regression Object
}
\description{
This class of objects is returned by the \code{smoothSurvReg} class of functions to represent
a fitted smoothed survival regression model.

Objects of this class have methods for the functions \code{print}, \code{summary}, \code{plot},
\code{residuals}, \code{survfit}.
}
\section{COMPONENTS I}{
The following components must be included in a legitimate \code{smoothSurvReg} object.
\describe{
\item{fail}{
Indicator of the failure of the fitting procedure. Possible values are
0 for no problems,
3 if the iteration process was stopped because of non-positive definite minus Hessian,
4 if the eiteration process was stopped because too many halving steps were performed,
5 if it was not possible to find the three reference knots (it was not then possible
to perform optimization with respect to the full parameter vector),
6 if the maximal number of iterations was performed without reaching a convergence.
The \code{fail} component is increased by 10 if the final minus Hessian of the penalized
log-likelihood was not positive definite.
The \code{fail} component is further increased by 20 if the computed effective 
degrees of freedom were non-positive.
The \code{fail} component is further increased by 40 if there are negative estimates
of standard errors for some regression parameters.
The \code{fail} component is 99 or higher if the fitting procedure failed at all and
there is no fit produced.
}
}
}
\section{COMPONENTS II}{
The following components must be included in a legitimate \code{smoothSurvReg} object
if \code{fail} is lower than 99.
\describe{
\item{regres}{
Estimates of the regression parameters \eqn{\alpha, \beta, \sigma}{alpha, beta, sigma} if
these have been estimated with their standard errors stored in a data frame
with \code{colnames} ``Value'', ``Std.Error'', ``Std.Error2'' and \code{rownames} derived from
the names of the design matrix with ``(Intercept)'' for the intercept, ``Scale'' for the scale
and ``Log(scale)'' for the log-scale. If the log-scale depends on covariates then rows
named ``LScale.(Intercept)'', ``LScale.cov1'' etc. give estimates of regression parameters
for log-scale.
The two standard errors are computed using either
\code{var} or \code{var2} described below.
}
\item{spline}{
Description of the fitted error density. 
A data frame with \code{colnames} ``Knot'', ``SD basis'',
``c coef.'', ``Std.Error.c'', ``Std.Error2.c'',
``a coef.'', ``Std.Error.a'' and ``Std.Error2.a''
and \code{rownames} knot[1], \dots, knot[g] where
\eqn{g}{g} stands for the number of basis G-splines. The column ``Knot'' contains the knots
in ascending order, ``SD basis'' the standard deviation of an appropriate basis G-spline,
``c coef.'' estimates of the G-spline coefficients and ``Std.Error.c'' and ``Std.Error2.c'' 
the estimates of their standard errors based either on \code{var} or \code{var2}. 
The column ``a coef.'' contains the estimates of transformed \eqn{c}{c} coefficients where 
\deqn{c_j = \frac{\exp(a_j)}{\sum_{l=1}^{g}\exp(a_l)}, j = 1,\dots, g.}{%
      c[j] = exp(a[j])/[exp(a[1]) + ... + exp(a[g])], j = 1, ..., g.}
If the error distribution is estimated, one of the \eqn{a}{a} coefficients is set to zero and  
two other \eqn{a}{a}'s are expressed as a function of the remaining \eqn{a}{a} coefficients
(to avoid equality constraints concerning the mean and the variance of the error distribution).
The standard error for these three \eqn{a}{a} coefficients is then not available (it is equal
to \code{NA}).
Standard error is set to
\code{NaN} is a diagonal element of the appropriate covariance matrix was negative.
}
\item{loglik}{
Maximized penalized log-likelihood, log-likelihood and the penalty term. A data frame with 
one row and three columns named ``Log Likelihood'', ``Penalty'' 
and ``Penalized Log Likelihood''.
}
\item{aic}{
Akaike's information criterion of the fitted model computed as a maximized value of the 
penalized log-likelihood minus the effective degrees of freedom.
}
\item{degree.smooth}{
Effective degrees of freedom, number of parameters and related information. A data frame 
with one row
and columns named ``Lambda'', ``Log(Lambda)'', ``df'', 
``Number of parameters'', ``Mean param.'', ``Scale param.'', ``Spline param.'' where
``Lambda'' gives the value of the tunning parameter used in the final (optimal) fit, 
``df'' the effective degrees of freedom, 
``Number of parameters'' the real number of parameters and 
``Mean param.'', ``Scale param.'' and ``Spline param.'' its decomposition. 
Note that if G-spline coefficients are estimated ``Spline param.''
is equal to the number of basis G-spline with non-zero coefficients minus three.
}
\item{var}{
The estimate of the covariance matrix of the estimates based on the Bayesian approximation.
It is equal to the inverse of the converged minus Hessian of the penalized log-likelihood.
Note that there are no columns and rows corresponding to the three
transformed G-spline coefficients since these are functions of the remaining transformed G-spline coefficients
(to avoid equality constraints).
}
\item{var2}{
The estimate of the covariance matrix of the estimates based on the asymptotic theory
for penalized models. It is equal to \eqn{H^{-1}\,I\,H^{-1}}{H^{-1} I H^{-1}} where
H is converged minus Hessian of the penalized log-likelihood and I is converged minus
Hessian of the log-likelihood component of the penalized log-likelihood.
}
\item{dCdD}{
A matrix with derivatives of \eqn{c}{c} spline coefficients with respect
to \eqn{d}{d} spline coefficients (these are \eqn{a}{a} coefficients
with three of them omitted). This matrix can be used later to compute
estimates and standard errors of functions of original parameters using
a Delta method. For closer definition of \eqn{d}{d} coefficients see an
enclosed document.
}  
\item{iter}{
Used number of iterations to fit the model with the optimal \eqn{\lambda}{lambda}.
}
\item{estimated}{
Indicator of what has really been estimated and not fixed. A four-component vector with 
component names ``(Intercept)'', ``Scale'', ``ccoef'', ``common.logscale''. 
The first component is TRUE if the intercept was included in the regression model.
The second component is TRUE if the scale parameter was not fixed, the third
component is TRUE is the G-spline coefficients were not fixed. The fourth
component is TRUE if the log-scale does not depend on covariates.
}
\item{warning}{
A data frame with one column called ``warnings'' and three rows called 
``Convergence'', ``Final minus Hessian'' and ``df'' 
containing a string information corresponding to the value of the \code{fail} component
of the object. It contains a string ``OK'' if there are no problems with the appropriate part
of the fitting process.
}
\item{H}{
Converged minus Hessian of the penalized log-likelihood. 
}
\item{I}{
Converged minus Hessian of the log-likelihood component of the penalized log-likelihood.
\eqn{I = H - G}{I = H - G}.
}
\item{G}{
Converged minus Hessian of the penalty term of the penalty term of the penalized log-likelihood.
\eqn{G = H - I}{G = H - I}.
}
\item{U}{
Converged score vector based on the penalized log-likelihood.
}
\item{na.action}{
The \code{na.action} attribute, if any, that was returned by the \code{na.action} routine.
}
\item{terms}{
The \code{terms} object used.
}
\item{formula}{
A symbolic description of the model to be fit.
}
\item{call}{
The matched call.
}
\item{init.dist}{
A string indicating the error distribution of the untransformed response to find the initial values. 
Possible values are ``lognormal'', ``loglogistic'', ``weibull''.
}
\item{model}{
If requested, the model frame used.
}
\item{x}{
The model matrix used.
}
\item{y}{
The response matrix used (two columns if there were no interval censored observations,
three columns if there were some interval censored observations). The last column
indicates the death status.
}
\item{z}{
The model matrix used for the expression of log-scale.
}
\item{init.spline}{
A data frame describing the initial error density. It has columns named ``Knot'', ``SD basis'',
``c coef.'' and rows named ``knot[1]'', \dots, ``knot[g]''.
}
\item{init.regres}{
Initial estimates of the regression parameters. A data frame with one column named ``Value''
and rows named as in the \code{regres} component of the \code{smoothSurvReg} object.
}
\item{adjust}{
Adjusted intercept and scale. A data frame with a column named ``Value'' and rows named ``(Intercept)''
and ``Scale''. ``(Intercept)'' gives the overall intercept taking into account the mean of the fitted
error distribution, ``Scale'' gives the overall scale taking into account the variance of the
fitted error distribution. If the error distribution is standardized
(always when G-spline coefficients are estimated)
then the ``(Intercept)'' is equal to the ``(Intercept)'' from the \code{regres} component and
``Scale'' is equal to the ``Scale'' of either \code{regres} or \code{init.regres} component.
\code{NA}'s appeare in this \code{data.frame} in the case that log-scale depends on covariates.
}
\item{error.dist}{
A data frame with columns named ``Mean'', ``Var'' and ``SD'' and a row named ``Error distribution:  ''
giving the mean, variance and the standard deviation of the fitted error distribution.
These are equal to 0, 1 and 1 if the G-spline coefficients were estimated.
}
\item{searched}{
Information concerning the searched values of the tunning paramater \eqn{\lambda}{lambda}
when looking for the best AIC. 
A data frame with columns named ``Lambda'', ``Log(Lambda)'', ``AIC'', ``df'',
``PenalLogLik'', ``LogLik'', ``nOfParm'', ``fail''.
}
}
}
\author{
\eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.ac.be}
}
\keyword{survival}
\keyword{smooth}
