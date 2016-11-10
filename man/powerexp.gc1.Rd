
\name{powerexp.gc1}
\alias{powerexp.gc1}
\title{The Powered Exponential Correlation Function of Class \code{\link{corr.gc1}}}
\usage{
powerexp.gc1(range = 0, kappa = 1, nugget = 0)
}
\arguments{
\item{range}{a non-negative scalar of the range parameter in powered
exponential correlation function.}

\item{kappa}{a scalar between 0 and 2; the value of the shape
parameter in the powered exponential correlation function.}

\item{nugget}{the nugget effect of the correlation function. A value between 0 and 1.}
}
\value{
An object of class \code{\link{corr.gc1}} representing the correlation component.
}
\description{
The powered exponential correlation function used for simulating geostatistical data
 from Gaussian copula models.


}
\details{
The powered exponential correlation function with a nugget \eqn{\tau^2} is of the form:

\deqn{\rho(h) = (1-\tau^2) exp((-h/\phi)^\kappa)}

  when \eqn{h > 0} and \eqn{\rho(h) = 1} when \eqn{h = 0}.
  Here \eqn{h} is distance, \eqn{\phi}{phi} is range parameter, \eqn{\kappa}{kappa} is the shape parameter and
  \eqn{\tau^2}{tau^2} is the nugget effect.

  When using the powered exponential correlation function, note that \eqn{0<\kappa \le 2}.
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{matern.gc1}},  \code{\link{spherical.gc1}}
}

\keyword{Correlation}% use one of  RShowDoc("KEYWORDS")
\keyword{Simulation}% __ONLY ONE__ keyword per line
