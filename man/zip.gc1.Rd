
\name{zip.gc1}
\alias{zip.gc1}
\title{The Zero-inflated Poisson Marginal of Class \code{\link{marginal.gc1}}}
\usage{
zip.gc1(mu = 1, od = 1)
}
\arguments{
\item{mu}{a non-negative scalar of the mean parameter.}

\item{od}{a non-negative scalar of the overdispersion parameter.}
}
\value{
An object of class \code{\link{marginal.gc1}} representing the marginal component.
}

\description{
The zero-inflated Poisson marginal used for simulation and computing
correlation in trans-Gaussian random field, parameterized in terms of the mean and overdispersion.
}

\details{
The zero-inflated Poisson distribution with parameters \code{mu = a} and \code{od = b} has density

\deqn{b/(1+b) + exp(-(a+ab))/(1+b)}
when \eqn{y = 0}, and

\deqn{exp(-(a+ab))*(a+ab)^y/((1+b)y!)}
when \eqn{y = 1, 2, \ldots}


Under this parameterization, \eqn{var(Y)= mu + od*mu^2}, where
\eqn{mu} is the mean parameter and \eqn{od} is the overdispersion parameter.
For more details see Han and De Oliveira (2016).
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}

\references{
Han, Z. and De Oliveira, V. (2016) On the correlation structure of Gaussian
            copula models for geostatistical count data.
            \emph{Australian and New Zealand Journal of Statistics, 58:47-69}.
}

\seealso{
\code{\link{marginal.gc1}},  \code{\link{beta.gc1}}, \code{\link{binomial.gc1}},
           \code{\link{gm.gc1}}, \code{\link{gaussian.gc1}},
           \code{\link{negbin.gc1}}, \code{\link{poisson.gc1}},
          \code{\link{weibull.gc1}}
}

\keyword{Marginal}
\keyword{Simulation}
