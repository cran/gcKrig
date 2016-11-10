
\name{gaussian.gc1}
\alias{gaussian.gc1}
\title{The Gaussian Marginal of Class \code{\link{marginal.gc1}}}
\usage{
gaussian.gc1(mean = 0, sd = 1)
}
\arguments{
\item{mean}{the mean of the Gaussian distribution, a scalar.}

\item{sd}{a positive scalar, the standard deviation of the Gaussian distribution.}
}
\value{
An object of class \code{\link{marginal.gc1}} representing the marginal component.
}
\description{
The Gaussian marginal used for simulation and computing correlation in
trans-Gaussian random field.
}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{marginal.gc1}}, \code{\link{beta.gc1}}, \code{\link{binomial.gc1}},
          \code{\link{gm.gc1}}, \code{\link{negbin.gc1}},
          \code{\link{poisson.gc1}}, \code{\link{weibull.gc1}}, \code{\link{zip.gc1}}
}

\keyword{Marginal}
\keyword{Simulation}
