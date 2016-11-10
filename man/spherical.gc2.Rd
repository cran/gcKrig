
\name{spherical.gc2}
\alias{spherical.gc2}
\title{The Spherical Correlation Function of Class \code{\link{corr.gc2}}}
\usage{
spherical.gc2(nugget = TRUE)
}
\arguments{
\item{nugget}{if \code{nugget = TRUE} then the correlation is estimated with a nugget effect between 0 and 1.
if \code{nugget = FALSE} then the correlation is estimated without the nugget effect.}
}

\value{
An object of class \code{\link{corr.gc2}} representing the correlation component.
}
\description{
The spherical correlation function used for likelihood inference and spatial prediction in Gaussian copula models.
}
\details{
See \code{\link{spherical.gc1}} for parameterization details.
}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{matern.gc2}},  \code{\link{powerexp.gc2}}
}

\keyword{Correlation}% use one of  RShowDoc("KEYWORDS")
\keyword{Estimation}% __ONLY ONE__ keyword per line
\keyword{Prediction}
