\name{marginal.gc2}
\alias{marginal.gc2}
\title{
Marginals for Likelihood Inference and Prediction in Gaussian Copula Models for Geostatistical Count Data
}
\description{
Class of marginals available in \code{gcKrig} library for likelihood inference and plug-in
prediction. In both cases the parameter values of the marginals are unknown thus need to be
estimated from the data.
}

\details{
This class is different from
the class \code{\link{marginal.gc1}}, in which the parameters are given by users.
Class \code{\link{marginal.gc2}} is used in the likelihood inference and spatial prediction, while
the class \code{\link{marginal.gc1}} is used in the simulation and correlation assessment.
The former case requires the input of the dataset while the latter requires the input of the parameters.
}
\value{
At the moment, the following marginals are implemented:
\tabular{ll}{

\code{\link{binomial.gc2}} \tab binomial marginals. \cr
\code{\link{negbin.gc2}} \tab negative binomial marginals. \cr
\code{\link{poisson.gc2}} \tab Poisson marginals. \cr
\code{\link{zip.gc2}} \tab zero-inflated Poisson marginals. \cr
}
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}

\seealso{\code{\link{binomial.gc2}},
         \code{\link{negbin.gc2}},
         \code{\link{poisson.gc2}},
         \code{\link{zip.gc2}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Marginal}
\keyword{Estimation}
\keyword{Prediction}
