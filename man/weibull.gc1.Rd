
\name{weibull.gc1}
\alias{weibull.gc1}
\title{The Weibull Marginal of Class \code{\link{marginal.gc1}}}
\usage{
weibull.gc1(shape = 1, scale = 1)
}
\arguments{
\item{shape}{a positive scalar of shape parameter in the Weibull distribution.}

\item{scale}{a positive scalar of scale parameter in the Weibull distribution.}
}
\value{
An object of class \code{\link{marginal.gc1}} representing the marginal component.
}
\description{
The Weibull marginal used for simulation and computing correlation in
trans-Gaussian random field.

The Weibull distribution with \code{shape} parameter \eqn{a} and \code{scale} parameter \eqn{b}
has density given by

\deqn{(a/b) (x/b)^{a-1} exp(-(x/b)^a) }{(a/b) (x/b)^(a-1) exp(-(x/b)^a)}
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{marginal.gc1}},  \code{\link{beta.gc1}},
          \code{\link{binomial.gc1}}, \code{\link{gm.gc1}},
          \code{\link{gaussian.gc1}}, \code{\link{negbin.gc1}},
          \code{\link{poisson.gc1}}, \code{\link{zip.gc1}}
}

\keyword{Marginal}
\keyword{Simulation}
