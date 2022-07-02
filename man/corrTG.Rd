
\name{corrTG}
\alias{corrTG}
\title{Compute the Correlation in  Transformed Gaussian Random Fields}
\usage{
corrTG(marg1, marg2, corrGauss = 0.5, method = "integral", nrep = 1000,
       kmax = 10, earlystop = FALSE, epscut = 1e-3)
}

\description{
	This function implements two general methods for computing the correlation function in
	a transformed Gaussian random field.
}


\arguments{
\item{marg1}{an object of class \code{\link{marginal.gc}} specifying the first marginal distribution.}

\item{marg2}{an object of class \code{\link{marginal.gc}} specifying the second marginal distribution.}

\item{corrGauss}{the correlation in the Gaussian random field. Should be a scalar between 0 and 1.}

\item{method}{the computation method of calculating correlation in the transformed Gaussian random field.
Can be either "integral" or "mc".
If use "integral" then a series expansion based on the Hermite Polynomials will be used to approximate the
correlation, see De Oliveira (2013) or Han and De Oliveira (2016).
If use "mc" then the Monte Carlo method will be used.}

\item{nrep}{the Monte Carlo size in computing the correlation. Only need to be specified if \code{method = "mc"}.}

\item{kmax}{the maximum number of terms used in the series summation (with Hermite polynomial expansion).
            Only need to be specified if \code{method = "integral"}.}

\item{earlystop}{whether or not to allow the series summation to stop early.
                 If \code{earlystop = FALSE} then a total number of kmax terms will be kept.
                 If \code{earlystop = TRUE} then the series will be automatically truncated if the absolute values
                 of the three consecutive terms are smaller than \code{epscut}.
                 }

\item{epscut}{a small positive value used to truncate the series.}

}
\value{
If \code{method = "mc"} the output is a scalar between 0 and 1, denoting the correlation of the transformed Gaussian random field. If \code{method = "integral"} the output is a scalar of the correlation in the transformed Gaussian random field, and a list of the values in the series expansion based on the integral with Hermite polynomials.
}
\examples{
\dontrun{
corrTG(marg1 = poisson.gc(lambda = 10), marg2 = binomial.gc(size = 1, prob = 0.1),
       corrGauss = 0.5, method = "integral", kmax = 10, earlystop = TRUE, epscut = 1e-5)
set.seed(12345)
corrTG(marg1 = poisson.gc(lambda = 10), marg2 = binomial.gc(size = 1, prob = 0.1),
       corrGauss = 0.5, nrep = 100000, method = "mc")
}
}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\references{
De Oliveira, V. (2013) Hierarchical Poisson models for spatial count data.
             \emph{Journal of Multivariate Analysis,122:393-408}.

Han, Z. and De Oliveira, V. (2016) On the correlation structure of Gaussian
            copula models for geostatistical count data.
            \emph{Australian and New Zealand Journal of Statistics, 58:47-69}.

Han, Z. and De Oliveira, V. (2018) gcKrig: An R Package for the Analysis of Geostatistical Count Data Using Gaussian Copulas.
\emph{Journal of Statistical Software}, \bold{87}(13), 1--32. \doi{10.18637/jss.v087.i13}.
}


\keyword{Correlation Structure}
