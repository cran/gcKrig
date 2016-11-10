\name{marginal.gc1}
\alias{marginal.gc1}
\title{
Marginals for Simulating Geostatistical Data from Gaussian Copula Models and
Correlation Matching for Trans-Gaussian Random Field (NORTA)
}
\description{
Class of marginals available in \code{gcKrig} library for simulation and
correlation structure assessment only. In both cases, the parameter values of
the marginals are given by users.
}

\details{
Parameter values need to be specified with the marginals.
For different marginals, parametrization may be different so the names of the
arguments may be different.


This class is different from
the class \code{\link{marginal.gc2}}, in which the parameters are unknown.
Class \code{\link{marginal.gc1}} is used in simulation and correlation assessment only.
Thus parameters need to be given by users.
Class \code{\link{marginal.gc2}} is used in likelihood inference and prediction thus parameters will be estimated.
}
\value{
At the moment, the following marginals are implemented:
\tabular{ll}{

\code{\link{beta.gc1}} \tab beta marginals. \cr
\code{\link{binomial.gc1}} \tab binomial marginals. \cr
\code{\link{gm.gc1}} \tab gamma marginals. \cr
\code{\link{gaussian.gc1}} \tab Gaussian marginals. \cr
\code{\link{negbin.gc1}} \tab negative binomial marginals. \cr
\code{\link{poisson.gc1}} \tab Poisson marginals. \cr
\code{\link{weibull.gc1}} \tab Weibull marginals. \cr
\code{\link{zip.gc1}} \tab zero-inflated Poisson marginals. \cr
}
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}


\seealso{\code{\link{beta.gc1}},
         \code{\link{binomial.gc1}},
         \code{\link{gm.gc1}},
         \code{\link{gaussian.gc1}},
         \code{\link{negbin.gc1}},
         \code{\link{poisson.gc1}},
         \code{\link{weibull.gc1}},
         \code{\link{zip.gc1}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Simulation}% use one of  RShowDoc("KEYWORDS")
\keyword{Marginal}
