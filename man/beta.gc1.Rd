
\name{beta.gc1}
\alias{beta.gc1}
\title{The Beta Marginal of Class \code{\link{marginal.gc1}}}
\usage{
beta.gc1(shape1 = 1, shape2 = 1)
}
\arguments{
\item{shape1}{non-negative scalar, the shape parameter of the Beta distribution.}

\item{shape2}{non-negative scalar, another shape parameter of the Beta distribution.}
}
\value{
An object of class \code{\link{marginal.gc1}} representing the marginal component.
}
\description{
The Beta marginal used for simulation and computing correlation in the
trans-Gaussian random field.
}

\details{
The Beta distribution with parameters \code{shape1 = a} and \code{shape2 = b} has density


\deqn{\frac{\Gamma(a+b)}{\Gamma(a) \Gamma(b)} \Big(\frac{b}{b+a}\Big)^b
        \Big(1 - \frac{b}{b+a}\Big)^y }{(\Gamma(a+b)/(\Gamma(a)\Gamma(b))) * x^(a-1) * (1-x)^(b-1)
}


for \eqn{a > 0, b > 0} and \eqn{0 \le x \le 1} where the boundary values at
\code{x = 0} or \code{x = 1} are defined as by continuity (as limits).
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{marginal.gc1}}, \code{\link{binomial.gc1}}, \code{\link{gm.gc1}},
          \code{\link{gaussian.gc1}}, \code{\link{negbin.gc1}},
          \code{\link{poisson.gc1}}, \code{\link{weibull.gc1}}, \code{\link{zip.gc1}}
}

\keyword{Marginal}
\keyword{Simulation}

