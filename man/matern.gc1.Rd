
\name{matern.gc1}
\alias{matern.gc1}
\title{The Matern Correlation Function of Class \code{\link{corr.gc1}}}
\usage{
matern.gc1(range = 0, kappa = 0.5, nugget = 0)
}
\arguments{
\item{range}{a non-negative scalar of the range parameter in Matern correlation function.}

\item{kappa}{a non-negative scalar of the shape parameter in the Matern correlation function.
The default kappa = 0.5 corresponds to an exponential correlation model.}

\item{nugget}{the nugget effect of the correlation function. A scalar between 0 and 1.}
}
\value{
An object of class \code{\link{corr.gc1}} representing the correlation component.
}
\description{
The Matern correlation function used for simulating geostatistical data
 from Gaussian copula models.
}
\details{
The Matern correlation function with a nugget \eqn{\tau^2} is of the form:

\deqn{\rho(h) =
    \frac{1-\tau^2}{2^{\kappa-1}\Gamma(\kappa)}\Big(\frac{h}{\phi}\Big)^\kappa
    K_{\kappa}\Big(\frac{h}{\phi}\Big)}{\rho(h) = (1-\tau^2)*
    (1/(2^(\kappa-1) * \Gamma(\kappa))) * ((h/\phi)^\kappa) * K_{\kappa}(h/\phi)}

  when \eqn{h > 0} and \eqn{\rho(h) = 1} when \eqn{h = 0}. Here
  \eqn{\phi}{phi} is range parameter, \eqn{\kappa}{kappa} is the shape parameter and
  \eqn{\tau^2}{tau^2} is the nugget parameter.
  \eqn{K_\kappa(\cdot)}{K_\kappa(\ldots)} denotes the modified Bessel function of the third
  kind of order \eqn{\kappa}{kappa}.
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\references{
Diggle, P. and Ribeiro, P.J. (2007) Model-based Geostatistics. \emph{Springer}.
}
\seealso{
\code{\link{powerexp.gc1}},  \code{\link{spherical.gc1}}
}

\keyword{Correlation}% use one of  RShowDoc("KEYWORDS")
\keyword{Simulation}% __ONLY ONE__ keyword per line
