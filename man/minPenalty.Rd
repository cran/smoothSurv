\name{minPenalty}
\alias{minPenalty}
\title{
  Minimize the penalty term under the two (mean and variance)
  constraints
}  
\description{
  This function minimizes \deqn{\frac{1}{2}\sum_{j=m+1}^g \Bigl(\Delta^m a_j\Bigr)^2}{%
    0.5*sum[j=m+1][g] (Delta^m a[j])^2} with respect to \eqn{a_1,\dots, a_g}{a[1],...,a[g]}
  under the constraints
  \deqn{\sum_{j=1}^g c_j \mu_j = 0}{sum[j=1][g] c[j]mu[j] = 0} and
  \deqn{\sum_{j=1}^g c_j (\mu_j^2 + \sigma_0^2) = 1,}{sum[j=1][g] c[j](mu[j]^2 + sigma0^2) = 1,}
  where 
  \deqn{c_j = \frac{\exp(a_j)}{\sum_{l=1}^g}\exp(a_l)}{c[j] = exp(a[j])/sum[l=1][g]exp(a[l])}
  with one of \eqn{a}{a}'s fixed to zero.

  Note that the minimum is always zero. We are thus mainly interested in
  the point where the minimum is reached.
}
\usage{
minPenalty(knots = NULL, dist.range = c(-6, 6), by.knots = 0.3, sdspline = NULL,
    difforder = 3, init.c,
    maxiter = 200, rel.tolerance = 1e-10, toler.chol = 1e-15, toler.eigen = 1e-3,
    maxhalf = 10, debug = 0, info = TRUE)
}
\arguments{
   \item{knots}{
       A vector of knots \eqn{\mu_1,\dots,\mu_g}{mu[1], ..., mu[g]}.
   }       
   \item{dist.range}{
       Approximate minimal and maximal knot. If not given by \code{knots} the knots
       are determined as \code{c(seq(0, dist.range[2], by = by.knots), seq(0, dist.range[1], by = -by.knots))}.
       The sequence of knots is sorted and multiple entries are removed.
   }
   \item{by.knots}{
       The distance between the two knots used when building a vector of knots if these
       are not given by \code{knots}.
   }
   \item{sdspline}{
       Standard deviation \eqn{\sigma_0^2}{sigma0^2} of the basis
       G-spline (here it appeares only in the variance constraint).
       If not given it is determined as 2/3 times the maximal distance between the two knots. If 
       \code{sdspline} >= 1 it is changed to 0.9 to be able to satisfy the constraints.
   }
   \item{difforder}{
       The order of the finite difference used in the penalty term.
   }
   \item{init.c}{
      Optional vector of the initial values for the G-spline
      coefficients c, all values must lie between 0 and 1
      and must sum up to 1.
   }
   \item{maxiter}{
       Maximum number of Newton-Raphson iterations.
   }
   \item{rel.tolerance}{
       (Relative) tolerance to declare the convergence. For this
       function, the convergence is declared if absolute value of the
       penalty is lower than \code{rel.tolerance} and if both
       constraints are satisfied up to \code{rel.tolerance}.
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
}
\value{
  A list with the components ``spline'', ``penalty'', ``warning'',
  ``fail''.
  \item{spline}{A data frame with columns named ``Knot'', ``SD basis'',
    ``c coef.'' and ``a coef.'' which gives the optimal values of
    \eqn{c_1,\dots, c_g}{c[1],...,c[g]} and 
    \eqn{a_1,\dots, a_g}{a[1],...,a[g]} in the latter two columns. This
    data.frame can be further worked out using the function \code{\link{eval.Gspline}}.}
  \item{penalty}{The value of the penalty term when declaring
    convergence.}
  \item{warning}{Possible warnings concerning the convergence.}
  \item{fail}{Failure indicator. It is zero if everything went OK.}
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\examples{
optimum <- minPenalty(knots=seq(-4.2, 4.2, by = 0.3), sdspline=0.2, difforder=3)
where <- optimum$spline
print(where)
show <- eval.Gspline(where, seq(-4.2, 4.2, by=0.05))
plot(show, type="l", bty="n", lwd=2)
}
\keyword{optimize}
