\name{derivative.cc3}
\alias{derivative.cc3}
\title{
   Work Function for 'smoothSurvReg', currently nowhere used
}
\description{
   Function to compute derivatives of all 'c' G-spline coefficients
   with respect to chosen (g - 3) coefficients such that
   the whole vector of g 'c' coefficients satisfies the constraints.
}
\usage{
derivative.cc3(knots, sdspline, last.three, all = TRUE)
}
\arguments{
   \item{knots}{
      A vector of G-spline knots \eqn{\mu}{mu}.
   }
   \item{sdspline}{
      Standard deviation \eqn{sigma_0}{sigma0} of the basis G-spline .
   }
   \item{last.three}{
      Indeces of the three 'c' G-spline coefficients which are
      expressed as a function of the remaining (g - 3)
      'c' G-spline coefficients such that the three constraints
      are satisfied.
      This must be a vector of length 3 with three different numbers from
      \code{1:length(knots)}.
   }
   \item{all}{
      If \code{TRUE}, matrix (g - 2) x g (there is one zero column) is returned.
      If \code{FALSE}, matrix (g - 2) x 3 is returned.
      The first row is always an intercept.
      See details.
   }
}
\details{
   To satisfy the three constraints
   \deqn{\sum_{j=1}^g c_j = 1,}{c[1] + ... c[g] = 1,}
   \deqn{\sum_{j=1}^g c_j \mu_j = 0,}{c[1]mu[1] + ... c[g]mu[g] = 0,}
   \deqn{\sum_{j=1}^g c_j \mu_j^2 = 1 - \sigma_0^2}{c[1]mu[1]^2 + ... c[g]mu[g]^2 = 1 - sigma0^2}
   imposed on the G-spline we can express the three 'c' coefficients as a function of
   the remaining \eqn{g - 3}{g - 3} 'c' coefficients in the following way.
   \deqn{c_{k} = \omega_{0,k} + \sum_{j\neq last.three}\omega_{j,k} c_j, %
         \qquad k \in last.three,}{%
         c[k] = omega[0,k] + sum[j != last.three]omega[j,k] c[j],   %
         k = last.three[1], last.three[2], last.three[3],}
    where \eqn{\omega}{omega} coefficients are a function of knots and G-spline
    standard deviation.
    If we denote \eqn{d}{d} the vector \code{c[-last.three]} this function computes
    derivatives of \eqn{c}{c} w.r.t. \eqn{d}{d} together
    with the intercept term used to compute \eqn{c}{c} from
    \eqn{d}{d}. This is actually a matrix of \eqn{\omega}{omega}
    coefficients. If we denote it as \eqn{\Omega}{Omega} then
    if \code{all == TRUE}
    \deqn{c = \Omega_{1,\cdot}^T +  \Omega_{-1,\cdot}^T d}{%
          c = t(Omega_[1,]) +  t(Omega[-1,]) d}
    and if \code{all == FALSE}
    \deqn{c[last.three] = \Omega_{1,\cdot}^T +  \Omega_{-1,\cdot}^T d.}{%
          c[last.three] = t(Omega_[1,]) +  t(Omega[-1,]) d.}
}
\value{
   A matrix with \eqn{\omega}{omega} coefficients.
}
\note{
WARNING: This function was primarily used inside \code{\link{smoothSurvReg}}.
Consequently, it has very few error checks on its input arguments.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.ac.be}
}
\keyword{utilities}
