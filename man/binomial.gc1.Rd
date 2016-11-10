
\name{binomial.gc1}
\alias{binomial.gc1}
\title{The Binomial Marginal of Class \code{\link{marginal.gc1}}}
\usage{
binomial.gc1(size = 1, prob = 0.5)
}
\arguments{
\item{size}{number of trials (zero or more).}

\item{prob}{probability of success on each trial.}
}
\value{
An object of class \code{\link{marginal.gc1}} representing the marginal component.
}
\description{
The binomial marginal used for simulation and computing correlation in the
trans-Gaussian random field.
}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{marginal.gc1}},  \code{\link{beta.gc1}},
           \code{\link{gm.gc1}}, \code{\link{gaussian.gc1}},
          \code{\link{negbin.gc1}}, \code{\link{poisson.gc1}},
          \code{\link{weibull.gc1}}, \code{\link{zip.gc1}}
}

\keyword{Marginal}
\keyword{Simulation}
