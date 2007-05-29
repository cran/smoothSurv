\name{smoothSurvReg.control}
\alias{smoothSurvReg.control}
\title{
   More Options for 'smoothSurvReg'
}
\description{
   This function checks and sets the fitting options for smoothSurvReg.
   Its arguments can be used instead of \dots in a call to smoothSurvReg.
}
\usage{
smoothSurvReg.control(est.c = TRUE, est.scale = TRUE,
   maxiter = 200, firstiter = 0, rel.tolerance = 5e-5, toler.chol = 1e-15, toler.eigen = 1e-3,
   maxhalf = 10, debug = 0, info = TRUE, lambda.use = 1.0, sdspline = NULL,
   difforder = 3, dist.range = c(-6, 6), by.knots = 0.3, knots = NULL, nsplines = NULL, last.three = NULL)
}
\arguments{
   \item{est.c}{
       If TRUE the G-spline coefficients are estimated. Otherwise, they are fixed
       to the values given by \code{init.c} parameter of \code{\link{smoothSurvReg}}.
   }
   \item{est.scale}{
       If TRUE the scale parameter \eqn{\sigma}{sigma} is estimated. Otherwise,
       it is fixed to the value given by \code{init.scale} parameter
       of \code{\link{smoothSurvReg}}.
   }
   \item{maxiter}{
       Maximum number of Newton-Raphson iterations.
   }
   \item{firstiter}{
       The index of the first iteration. This option comes from older versions
       of this function.
   }
   \item{rel.tolerance}{
       (Relative) tolerance to declare the convergence. In this version of the function,
       the convergence is declared if the relative difference between two consecutive values
       of the penalized log-likelihood are smaller than rel.tolerance.
   }
   \item{toler.chol}{
       Tolerance to declare Cholesky decomposition singular.
   }
   \item{toler.eigen}{
       Tolerance to declare an eigen value of a matrix to be zero.
   }
   \item{maxhalf}{
       Maximum number of step-halving steps if updated estimate leads to a decrease
       of the objective function.
   }
   \item{debug}{
       If non-zero print debugging information.
   }
   \item{info}{
       If TRUE information concerning the iteration process is printed
       during the computation to the standard output.
   }
   \item{lambda.use}{
       The value of the tuning (penalty) parameter \eqn{\lambda}{\lambda} used 
       in a current fit by the \code{smoothSurvReg.fit} function.
       Value of this option is not interesting for the user. The parameter
       \code{lambda} of the function \code{smoothSurvReg} is more important for 
       the user.
   }
   \item{sdspline}{
       Standard deviation of the basis G-spline. If not given it is determined
       as 2/3 times the maximal distance between the two knots. If \code{est.c} = TRUE
       and \code{sdspline} >= 1 it is changed to 0.9 to be able to satisfy the constraints
       imposed to the fitted error distribution.
   }
   \item{difforder}{
       The order of the finite difference used in the penalty term.
   }
   \item{dist.range}{
       Approximate minimal and maximal knot. If not given by \code{knots} the knots
       are determined as \code{c(seq(0, dist.range[2], by = by.knots), seq(0, dist.range[1], by = -by.knots))}.
       The sequence of knots is sorted and multiple entries are removed.
   }
   \item{by.knots}{
       The distance between the two knots used when building a vector of knots if these
       are not given by \code{knots}. This option is ignored if \code{nsplines} is not \code{NULL}.
   }
   \item{knots}{
       A vector of knots.
   }
   \item{nsplines}{
       This option is ignored at this moment. It is used to give the number of G-splines
       to the function \code{\link{smoothSurvReg.fit}}.
   }
   \item{last.three}{
       A vector of length 3 with indeces of reference knots. The 'a' coefficient of
       the \code{knot[last.three[1]]} is then equal to zero, 'a' coefficients
       with indeces \code{last.three[2:3]} are expressed as a function of remaining
       'a' coefficients such that resulting error distribution has zero mean and unit
       variance. If \code{maxiter > 0} \code{last.three} is determined after
       the convergence is reached. If \code{maxiter == 0} \code{last.three} is used 
       to compute variance matrices.
   }
}
\value{
   A list with the same elements as the input except \code{dist.range} and \code{by.knots}
   is returned.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek@med.kuleuven.be}
}
\keyword{survival}
\keyword{smooth}
