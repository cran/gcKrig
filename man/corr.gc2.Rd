\name{corr.gc2}
\alias{corr.gc2}
\title{
Spatial Correlation Functions for Likelihood Inference and Spatial Prediction in Gaussian Copula Models
}
\description{
Class of isotropic correlation functions available in \code{gcKrig} library for likelihood
inference and spatial prediction.
}

\details{
By default, the nugget effect is included. To specify the model without the nugget effect
use \code{nugget = FALSE}.

For different correlation functions, parametrization may be different.

This class is different from
the class \code{\link{corr.gc1}}, in which the parameters are given by users.
Class \code{\link{corr.gc2}} is used for the purpose of estimation and prediction,
therefore, parameters are estimated from dataset.
}


\value{
At the moment, the following three correlation functins are implemented:
\tabular{ll}{

\code{\link{matern.gc2}} \tab the Matern correlation function. \cr
\code{\link{powerexp.gc2}} \tab the powered exponential correlation function. \cr
\code{\link{spherical.gc2}} \tab the spherical correlation function. \cr
}
}


\author{
Zifei Han \email{hanzifei1@gmail.com}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
         \code{\link{matern.gc2}},
         \code{\link{powerexp.gc2}},
         \code{\link{spherical.gc2}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Correlation Function}% use one of  RShowDoc("KEYWORDS")
\keyword{Estimation}% __ONLY ONE__ keyword per line
\keyword{Prediction}
