\name{derivative.expAD}
\alias{derivative.expAD}
\title{
   Work Function for 'smoothSurvReg', currently nowhere used
}
\description{
   Function to compute derivatives of exp(a) w.r.t. exp(d)
   where d stands for a shorter vector of 'a' G-spline coefficients.
}
\usage{
derivative.expAD(knots, sdspline, last.three, all = TRUE)
}
\arguments{
   \item{knots}{
      A vector of G-spline knots \eqn{\mu}{mu}.
   }
   \item{sdspline}{
      Standard deviation \eqn{\sigma_0}{sigma0} of the basis G-spline.
   }
   \item{last.three}{
      Indeces of the three 'a' G-spline coefficients which are
      expressed as a function of the remaining (g-3)
      'a' G-spline coefficients such that the three constraints
      are satisfied.
      This must be a vector of length 3 with three different numbers from
      \code{1:length(mu)}.
      Consequently, \code{a[last.three[1]] <- 0}.
   }
   \item{all}{
      If \code{TRUE}, matrix (g - 2) x g (there is one zero column) is returned.
      If \code{FALSE}, matrix (g - 2) x 2 is returned.
      The first row is always an intercept.
      See details.
   }
}
\details{
   To satisfy the three constraints
   \deqn{\sum_{j=1}^g c_j = 1,}{c[1] + ... c[g] = 1,}
   \deqn{\sum_{j=1}^g c_j \mu_j = 0,}{c[1]mu[1] + ... c[g]mu[g] = 0,}
   \deqn{\sum_{j=1}^g c_j \mu_j^2 = 1 - \sigma_0^2}{c[1]mu[1]^2 + ... c[g]mu[g]^2 = 1 - sigma0^2}
   imposed on the G-spline we use the following parametrization:
   \deqn{c_j = \frac{\exp(a_j)}{\sum_{l=1}^{g}\exp(a_l)}, j = 1,\dots, g.}{%
         c[j] = exp(a[j])/[exp(a[1]) + ... + exp(a[g])], j = 1, ..., g.}
   The constraints can be solved such that \code{a[last.three[1]] = 0} and
   \code{a[last.three[2:3]]} are expressed as a function of
   \code{a[-last.three]} in the following way:
   \deqn{a_{k} = \log\Bigl\{\omega_{0,k} + \sum_{j\neq last.three}\omega_{j,k}\exp(a_j)\Bigr\},%
         \qquad k = last.three[2], last.three[3],}{%
         a[k] = log[omega[0,k] + sum[j != last.three]omega[j,k]exp(a[j])],   %
         k = last.three[2], last.three[3],}
    where \eqn{\omega}{omega} coefficients are a function of knots and G-spline
    standard deviation.
    If we denote \eqn{d}{d} the vector \code{a[-last.three]} this function computes
    derivatives of \eqn{\exp(a)}{exp(a)} w.r.t. \eqn{\exp(d)}{exp(d)} together
    with the intercept term used to compute \eqn{\exp(a)}{exp(a)} from
    \eqn{\exp(d)}{exp(d)}. This is actually a matrix of \eqn{\omega}{omega}
    coefficients. If we denote it as \eqn{\Omega}{Omega} then
    if \code{all == TRUE}
    \deqn{\exp(a) = \Omega_{1,\cdot}^T +  \Omega_{-1,\cdot}^T\exp(d)}{%
          exp(a) = t(Omega_[1,]) +  t(Omega[-1,])exp(d)}
    and if \code{all == FALSE}
    \deqn{\exp(a[last.three[2:3]]) = \Omega_{1,\cdot}^T +  \Omega_{-1,\cdot}^T\exp(d).}{%
          exp(a[last.three[2:3]]) = t(Omega_[1,]) +  t(Omega[-1,])exp(d).}
}
\value{
   A matrix with \eqn{\omega}{omega} coefficients.
}
\note{
WARNING: This function is primarily used inside \code{\link{smoothSurvReg}}.
Consequently, it has very few error checks on its input arguments.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{komarek@karlin.mff.cuni.cz}
}
\keyword{internal}
\keyword{utilities}
