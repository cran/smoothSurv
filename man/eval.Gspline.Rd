\name{eval.Gspline}
\alias{eval.Gspline}
\title{
  Evaluate a G-spline in a grid of values
}
\description{
  This function computes values of
  \deqn{f(x) = \sum_{j=1}^g c_j \varphi_{\mu_j, \sigma_j^2}(x)}{%
    f(x) = sum[j=1][g] c[j]*phi(x; mu[j],sigma[j]^2)}
  in a grid of \eqn{x}{x} values.

  In above expression, \eqn{\varphi_{\mu_j, \sigma_j^2}(x)}{phi(x; mu[j],sigma[j]^2)}  
  denotes a density of \eqn{N(\mu_j, \sigma_j^2)}{N(mu[j], sigma[j]^2)}.
}  
\usage{
  eval.Gspline(Gspline, grid)
}
\arguments{
  \item{Gspline}{A data frame with at least three columns named
    ``Knot'', ``SD basis'' and ``c coef.'' which determine
    \eqn{\mu_1, \dots,\mu_g}{mu[1], ..., mu[g]},
    \eqn{\sigma_1, \dots, \sigma_g}{sigma[1], ..., sigma[g]} and
    \eqn{c_1,\dots, c_g}{c[1], ..., c[g]}. Data.frame with such
    properties can be found e.g. as \code{spline} component of the
    resulting object returned by functions \code{\link{smoothSurvReg}}
    and \code{\link{minPenalty}}.    
  }
  \item{grid}{A numeric vector giving the grid of \eqn{x}{x} values at
    which the G-spline is to be evaluated.
  }      
}
\value{A data.frame with columns named ``x'' (grid)  and ``y'' (G-spline
  values).
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}
\examples{
  spline <- minPenalty(knots=seq(-4.2, 4.2, by=0.3), sdspline=0.2, difforder=3)$spline
  values <- eval.Gspline(spline, seq(-4.5, 4.5, by=0.05))
  plot(values, type="l", bty="n", lwd=3)
}
\keyword{dplot}
