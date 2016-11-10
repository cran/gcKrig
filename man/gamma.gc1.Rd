
\name{gm.gc1}
\alias{gm.gc1}
\title{The Gamma Marginal of Class \code{\link{marginal.gc1}}}
\usage{
gm.gc1(shape = 1, rate = 1)
}
\arguments{
\item{shape}{a non-negative scalar, shape parameter of the Gamma distribution.}

\item{rate}{a non-negative scalar, rate parameter of the Gamma distribution.}
}
\value{
An object of class \code{\link{marginal.gc1}} representing the marginal component.
}
\description{
The Gamma marginal used for simulation and computing correlation in
trans-Gaussian random field.
}
\details{
The Gamma distribution with parameters \code{shape = a} and \code{rate = r} has density

\deqn{\frac{r^{a}}{\Gamma(a)} x^{a-1} exp(-rx)}{r^a/\Gamma(a) x^(a-1) exp(-rx)}
for \eqn{x \ge 0}, \eqn{a > 0} and \eqn{s > 0}.
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{marginal.gc1}}, \code{\link{beta.gc1}}, \code{\link{binomial.gc1}},
          \code{\link{gaussian.gc1}}, \code{\link{negbin.gc1}},
          \code{\link{poisson.gc1}}, \code{\link{weibull.gc1}}, \code{\link{zip.gc1}}
}

\keyword{Marginal}
\keyword{Simulation}
