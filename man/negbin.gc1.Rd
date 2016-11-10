
\name{negbin.gc1}
\alias{negbin.gc1}
\title{The Negative Binomial Marginal of Class \code{\link{marginal.gc1}}}
\usage{
negbin.gc1(mu = 1, od = 1)
}
\arguments{
\item{mu}{a non-negative scalar of the mean parameter.}

\item{od}{a non-negative scalar of the overdispersion parameter.}
}
\value{
An object of class \code{\link{marginal.gc1}} representing the marginal component.
}
\description{
The negative binomial marginal used for simulation and computing correlation in
trans-Gaussian random field, parameterized in terms of the mean and overdispersion.
}

\details{
The negative binomial distribution with parameters \code{mu = a} and \code{od = 1/b} has density

\deqn{\frac{\Gamma(y+b)}{\Gamma(b)y!} \Big(\frac{b}{b+a}\Big)^b
        \Big(1 - \frac{b}{b+a}\Big)^y }{\Gamma(y+b)/(\Gamma(b)y!) * (b/(b+a))^b * (a/(b+a))^y
}
which is called NB2 by Cameron and Trivedi (2013).
Under this parameterization, \eqn{var(Y)= mu + od*mu^2}, where
\eqn{mu} is the mean parameter and \eqn{od} is the overdispersion parameter.
For more details see Han and De Oliveira (2016).
}


\author{
Zifei Han \email{hanzifei1@gmail.com}
}

\references{
 Cameron,A.C. and Trivedi,P.K. (2013) Regression Analysis of Count Data.
 \emph{Cambridge University Press, 2nd Edition}.

 Han, Z. and De Oliveira, V. (2016) On the correlation structure of Gaussian
            copula models for geostatistical count data.
            \emph{Australian and New Zealand Journal of Statistics, 58:47-69}.
}


\seealso{
\code{\link{marginal.gc1}},  \code{\link{beta.gc1}},
          \code{\link{binomial.gc1}}, \code{\link{gm.gc1}},
          \code{\link{gaussian.gc1}}, \code{\link{poisson.gc1}},
          \code{\link{weibull.gc1}}, \code{\link{zip.gc1}}
}

\keyword{Marginal}
\keyword{Simulation}
