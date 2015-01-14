\name{beed.confband}
\alias{beed.confband}
\title{Nonparametric Bootstrap Confidence Intervals}
\description{
  Compute nonparametric bootstrap \eqn{(1-\alpha)\%} confidence bands for the Pickands dependence function.
}
\usage{
  beed.confband(data, x, d = 3, est = c('ht','md','cfg'), 
  margin = c('emp','Gev'), k = 13, nboot = 200, y=NULL, conf = 0.95, 
  matrix = FALSE, plot = FALSE, print = FALSE)
}
\arguments{
  \item{data}{A \eqn{(n \times d)}{(n x d)} matrix of component-wise maxima.}
  \item{x}{A \eqn{(m \times d)}{(m x d)} design matrix (see \bold{Details}).}
  \item{d}{A postive integer (greater than or equal to two) indicating the number of variables. 
  The trivariate case \code{d = 3} is the default.}
  \item{est}{A string denoting the estimation method (see \bold{Details}).}
    \item{margin}{A string denoting the type marginal distributions (see \bold{Details}).}
  \item{k}{A postive integer denoting the order of the Bernstein polynomial. \code{k = 13} is set by default.}
  \item{nboot}{A postive integer indicating the number of bootstrap replicates.}
  \item{y}{A numeric vector (of size \code{m}) with an initial estimate of the Pickands function. If \code{NULL}, The initial estimation is performed by using the estimation method chosen in \code{est}.}
  \item{conf}{A real value in \eqn{(0,1)} denoting the confidence level of the interval. The value \code{conf = 0.95} is the default.}
  \item{matrix}{Logical; \code{FALSE} by default. If \code{TRUE}, and the dimension \code{d} is three (the default dimension), the value of \eqn{A} is collected in a square matrix.}
  \item{plot}{Logical; \code{FALSE} by default. If \code{TRUE}, the confidence bands are plotted.}
  \item{print}{Logical; \code{FALSE} by default. If \code{TRUE}, the number of the iteration is printed.}
}

\details{
Two methods for computing bootstrap \eqn{(1-\alpha)\%} point-wise and simultaneous confident bands for the Pickands dependence function are used. 

The first method derives the confidence bands computing the point-wise \eqn{\alpha/2}
and \eqn{1-\alpha/2} quantiles of the bootstrap sample distribution of the Pickands dependence Bernstein based estimator.

The second method derives the confidence bands, first  computing the point-wise \eqn{\alpha/2}
and \eqn{1-\alpha/2} quantiles of the bootstrap sample distribution of polynomial coefficients estimators, and then
the Pickands dependence is computed using the Bernstein polynomial representation. See Marcon et al. (2014) for details.

 Most of the settings are the same as in the function \code{\link{beed}}.
}

\note{
This routine relies on the bootstrap routine (see \code{\link{beed.boot}}).
}
  
  \value{
    \item{A}{Estimate of the Pickands dependence function.}
    \item{bootA}{A matrix with \code{nboot} columns that reports the estimates of the Pickands function for each data resampling.}
    \item{A.up.beta/A.low.beta}{Vectors of upper and lower bands of the Pickands dependence function obtained using the bootstrap sampling distribution of the polynomial coefficients estimator.}
    \item{A.up.pointwise/A.low.pointwise}{Vectors of upper and lower bands of the Pickands dependence function obtained using the bootstrap sampling distribution of the Pickands dependence function estimator.}
    \item{up.beta/low.beta}{Vectors of upper and lower bounds of the bootstrap sampling distribution of the polynomial coefficients estimator.}
  }
  
  \references{
    Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. (2014)
    Multivariate Nonparametric Estimation of the Pickands   
    Dependence Function using Bernstein Polynomials.
    \emph{Submitted}.
  }
  
  \seealso{\code{\link{beed}}, \code{\link{beed.boot}}.}
  

\examples{
w <- seq(0, 1, length = 100)  
data <- rbvevd(50, dep = 0.4, model = 'log', mar1 = c(1,1,1))

# Note you should consider 500 bootstrap replications. 
# In order to obtain fastest results we used 30!
cb <- beed.confband(data, cbind(w, 1-w), 2, 'md', 'emp', 20, 30) 
                    
plot(w, w, type='n', xlab="w", ylab="A(w)", ylim=c(.5,1))
polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty=1, lwd=1, border='grey')
lines(w, cb$A, lty=1)
lines(w, cb$A.up.beta, lty=2)
lines(w, cb$A.low.beta, lty=2)
}
  
\keyword{Bootstrap confidence bands}