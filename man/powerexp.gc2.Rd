
\name{powerexp.gc2}
\alias{powerexp.gc2}
\title{The Powered Exponential Correlation Function of Class \code{\link{corr.gc2}}}
\usage{
powerexp.gc2(kappa = 1, nugget = TRUE)
}
\arguments{
\item{kappa}{a scalar between 0 and 2; the value of the shape
parameter in the powered exponential correlation function.}

\item{nugget}{if \code{nugget = TRUE} then the correlation is estimated with a nugget effect between 0 and 1.
if \code{nugget = FALSE} then the correlation is estimated without the nugget effect.}


}
\value{
An object of class \code{\link{corr.gc2}} representing the correlation component.
}
\description{
The powered exponential correlation function used for likelihood inference and spatial prediction in
Gaussian copula models.
}
\details{
See \code{\link{powerexp.gc1}} for parameterization details.

Note that parameter \code{kappa} is given by the users, not estimated by the program.
}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{matern.gc2}},  \code{\link{spherical.gc2}}
}

\keyword{Correlation}% use one of  RShowDoc("KEYWORDS")
\keyword{Estimation}% __ONLY ONE__ keyword per line
\keyword{Prediction}
