\name{negbin.gc2}
\alias{binomial.gc2}
\alias{negbin.gc2}
\alias{poisson.gc2}
\alias{zip.gc2}

\title{
Marginals for Likelihood Inference and Prediction in Gaussian Copula Models for Geostatistical Count Data
}
\description{
These functions set the marginals in Gaussian copula models for likelihood inference and
spatial prediction.
}
\usage{
binomial.gc2(link = "logit")
negbin.gc2(link = "log")
poisson.gc2(link = "log")
zip.gc2(link = "log")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{the model link function.}
}

\value{
An object of class \code{\link{marginal.gc2}} representing the marginal component of the model.
For parameterization details see \code{\link{negbin.gc1}} and \code{\link{zip.gc1}}.
}
\references{
Cribari-Neto, F. and Zeileis, A. (2010) Beta regression in R. \emph{Journal of Statistical Software, 34:1-24}. \url{http://www.jstatsoft.org/v34/i02/}.

Ferrari, S.L.P. and Cribari-Neto, F. (2004) Beta regression for modeling rates and proportions. \emph{Journal of Applied Statistics, 31:799-815}.

Masarotto, G. and Varin, C. (2012) Gaussian copula marginal regression.
\emph{Electronic Journal of Statistics, 6:1517-1549}.
\url{https://projecteuclid.org/euclid.ejs/1346421603}.

Masarotto, G. and Varin, C. (2016) Gaussian Copula Regression in R. \emph{Journal of Statistical Software, To appear}.
\url{http://cristianovarin.weebly.com/uploads/1/5/1/5/15156956/masarotto_varin_gcmr-vignette.pdf}.
}


\author{
Zifei Han \email{hanzifei1@gmail.com}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Marginal}
\keyword{Estimation}
\keyword{Prediction}
