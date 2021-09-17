\name{beed.confband}
\alias{beed.confband}
\title{Nonparametric Bootstrap Confidence Intervals}
\description{
  Computes nonparametric bootstrap \eqn{(1-\alpha)\%} confidence bands for the Pickands dependence function.
}
\usage{
  beed.confband(data, x, d = 3, est = c("ht","md","cfg","pick"),
     margin=c("emp", "est", "exp", "frechet", "gumbel"), k = 13,
     nboot = 500, y = NULL, conf = 0.95, plot = FALSE, print = FALSE)
}
\arguments{
  \item{data}{\eqn{(n \times d)}{(n x d)} matrix of component-wise maxima.}
  \item{x}{\eqn{(m \times d)}{(m x d)} design matrix (see \bold{Details}).}
  \item{d}{postive integer (greater than or equal to two) indicating the number of variables (\code{d=3} by default).}
  \item{est}{string denoting the estimation method (see \bold{Details}).}
    \item{margin}{string denoting the type marginal distributions (see \bold{Details}).}
  \item{k}{postive integer denoting the order of the Bernstein polynomial (\code{k=13} by default).}
  \item{nboot}{postive integer indicating the number of bootstrap replicates.}
  \item{y}{numeric vector (of size \code{m}) with an initial estimate of the Pickands function. If \code{NULL}, the initial estimation is performed by using the estimation method chosen in \code{est}.}
  \item{conf}{real value in \eqn{(0,1)} denoting the confidence level of the interval. The value \code{conf=0.95} is the default.}
  \item{plot}{logical; \code{FALSE} by default. If \code{TRUE}, the confidence bands are plotted.}
  \item{print}{logical; \code{FALSE} by default. If \code{TRUE}, the number of the iteration is printed.}
}

\details{
Two methods for computing bootstrap \eqn{(1-\alpha)\%} point-wise and simultaneous confidence bands for the Pickands dependence function are used.

The first method derives the confidence bands computing the point-wise \eqn{\alpha/2} and \eqn{1-\alpha/2} quantiles of the bootstrap sample distribution of the Pickands dependence Bernstein based estimator.

The second method derives the confidence bands, first  computing the point-wise \eqn{\alpha/2} and \eqn{1-\alpha/2} quantiles of the bootstrap sample distribution of polynomial coefficient estimators, and then the Pickands dependence is computed using the Bernstein polynomial representation. See Marcon et al. (2017) for details.

Most of the settings are the same as in the function \code{\link{beed}}.
}

\note{
This routine relies on the bootstrap routine (see \code{\link{beed.boot}}).
}

  \value{
    \item{A}{numeric vector of the Pickands dependence function estimated.}
    \item{bootA}{matrix with \code{nboot} columns that reports the estimates of the Pickands function for each data resampling.}
    \item{A.up.beta/A.low.beta}{vectors of upper and lower bands of the Pickands dependence function obtained using the bootstrap sampling distribution of the polynomial coefficients estimator.}
    \item{A.up.pointwise/A.low.pointwise}{vectors of upper and lower bands of the Pickands dependence function obtained using the bootstrap sampling distribution of the Pickands dependence function estimator.}
    \item{up.beta/low.beta}{vectors of upper and lower bounds of the bootstrap sampling distribution of the polynomial coefficients estimator.}
  }

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://mypage.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com/};
	Giulia Marcon, \email{giuliamarcongm@gmail.com}
}


\references{
  Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. and Segers, J. (2017)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials.
  \emph{Journal of Statistical Planning and Inference}, 183, 1-17.
}

\seealso{\code{\link{beed}}, \code{\link{beed.boot}}.}


\examples{
if (interactive()){

x <- ExtremalDep:::simplex(2)
data <- rbvevd(50, dep = 0.4, model = "log", mar1 = c(1,1,1))

# Note you should consider 500 bootstrap replications.
# In order to obtain fastest results we used 50!
cb <- beed.confband(data, x, 2, "md", "emp", 20, 50, plot=TRUE)
}
}

\keyword{Bootstrap confidence bands}
