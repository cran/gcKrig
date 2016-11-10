
\name{corrTG}
\alias{corrTG}
\title{Compute the Correlation in  Transformed Gaussian Random Fields}
\usage{
corrTG(marg1, marg2, corrGaus = 0.5, method = "integral", nrep = 1000)
}

\description{
	This function implements two general methods for computing the correlation function in 
	a transformed Gaussian random field. 
}


\arguments{
\item{marg1}{an object of class \code{\link{marginal.gc1}} specifying the first marginal distribution.}

\item{marg2}{an object of class \code{\link{marginal.gc1}} specifying the second marginal distribution.}

\item{corrGaus}{the correlation in the Gaussian random field. Should be a scalar between 0 and 1.}

\item{method}{the computation method of calculating correlation in the transformed Gaussian random field.
Can be either "integral" or "mc".
If use "integral" then a series expansion based on the Hermite Polynomials will be used to approximate the
correlation, see De Oliveira (2013).
If use "mc" then the Monte Carlo method will be used.}

\item{nrep}{the Monte Carlo size in computing the correlation. Only need to be specified if \code{method = "mc"}.}
}
\value{
A scalar between 0 and 1, denoting the correlation of the transformed Gaussian random field.
}
\examples{
corrTG(marg1 = poisson.gc1(lambda = 10), marg2 = binomial.gc1(size = 1, prob = 0.1),
       corrGaus = 0.5, method = "integral")
set.seed(12345)
corrTG(marg1 = poisson.gc1(lambda = 10), marg2 = binomial.gc1(size = 1, prob = 0.1),
       corrGaus = 0.5, nrep = 100000, method = "mc")

}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\references{
De Oliveira, V. (2013) Hierarchical Poisson models for spatial count data.
             \emph{Journal of Multivariate Analysis,122:393-408}.
}


\keyword{Correlation Structure}
