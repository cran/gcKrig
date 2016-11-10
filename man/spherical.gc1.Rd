
\name{spherical.gc1}
\alias{spherical.gc1}
\title{The Spherical Correlation Function of Class \code{\link{corr.gc1}}}
\usage{
spherical.gc1(range = 0, nugget = 0)
}
\arguments{
\item{range}{a non-negative scalar of the range parameter in the spherical
correlation function.}

\item{nugget}{the nugget effect of the correlation function. A value between 0 and 1.}
}
\value{
An object of class \code{\link{corr.gc1}} representing the correlation component.
}
\description{
The spherical correlation function used for simulating geostatistical data
 from Gaussian copula models.
}
\details{
The spherical correlation function with a nugget \eqn{\tau^2} is of the form:

\deqn{\rho(h) = (1-\tau^2) (1 - 1.5(h/\phi) + 0.5(-h/\phi)^3) }
 when \eqn{h > 0} and \eqn{\rho(h) = 1} when \eqn{h = 0}, \eqn{h} is distance.
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{matern.gc1}},  \code{\link{spherical.gc1}}
}

\keyword{Correlation}% use one of  RShowDoc("KEYWORDS")
\keyword{Simulation}% __ONLY ONE__ keyword per line
