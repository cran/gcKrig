\name{corr.gc1}
\alias{corr.gc1}
\title{
Spatial Correlation Functions for Simulating Geostatistical Data from Gaussian Copula Models
}
\description{
Class of isotropic correlation functions available in the \code{gcKrig} library for simulation only.
The parameter values in the correlation functions should be specified by users.
}

\details{
Parameter values need to be specified with the correlation functions of this class.
For different correlation functions, parametrization may be different so the names of the
arguments may be different.

This class is different from
the class \code{\link{corr.gc2}}, in which the parameters are unknown.
Class \code{\link{corr.gc1}} is used for the purpose of simulation only.
Class \code{\link{corr.gc2}} is used in likelihood inference and prediction, therefore,
parameters will be estimated.
}


\value{
At the moment, the following three correlation functins are implemented:
\tabular{ll}{

\code{\link{matern.gc1}} \tab the Matern correlation function. \cr
\code{\link{powerexp.gc1}} \tab the powered exponential correlation function. \cr
\code{\link{spherical.gc1}} \tab the spherical correlation function. \cr
}
}


\author{
Zifei Han \email{hanzifei1@gmail.com}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
         \code{\link{matern.gc1}},
         \code{\link{powerexp.gc1}},
         \code{\link{spherical.gc1}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Correlation}% use one of  RShowDoc("KEYWORDS")
\keyword{Simulation}% __ONLY ONE__ keyword per line
