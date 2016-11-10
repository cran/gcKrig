
\name{matern.gc2}
\alias{matern.gc2}
\title{The Matern Correlation Function of Class \code{\link{corr.gc2}}}
\usage{
matern.gc2(kappa = 0.5, nugget = TRUE)
}
\arguments{
\item{kappa}{a non-negative scalar of the shape parameter in the Matern correlation function.
The default kappa = 0.5 corresponds to an exponential correlation model.}

\item{nugget}{if \code{nugget = TRUE} then the correlation is estimated with a nugget effect between 0 and 1.
if \code{nugget = FALSE} then the correlation is estimated without the nugget effect.}


}
\value{
An object of class \code{\link{corr.gc2}} representing the correlation component.
}
\description{
The Matern correlation function used for likelihood inference and spatial prediction in Gaussian copula models.
}
\details{
See \code{\link{matern.gc1}} for parameterization details.

Note that parameter \code{kappa} is given by the users, not estimated by the program.
}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\seealso{
\code{\link{powerexp.gc2}},  \code{\link{spherical.gc2}}
}

\keyword{Correlation}% use one of  RShowDoc("KEYWORDS")
\keyword{Estimation}% __ONLY ONE__ keyword per line
\keyword{Prediction}
