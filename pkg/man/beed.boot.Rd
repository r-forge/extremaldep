\name{beed.boot}
\alias{beed.boot}
\title{Bootstrap Resampling and Bernstein Estimation of Extremal Dependence}
\description{
 Compute \code{nboot} estimates of the uni- or multi-variate Pickands dependence function (using the 
 Bernstein polynomials approximation method) on the basis of the bootstrap resampling of the data.
}
\usage{
  beed.boot(data, x, d = 3, est = c('ht','md','cfg'), margin = c('emp','Gev'), 
            k = 13, nboot = 200, y = NULL, matrix = FALSE, print = FALSE)
}
\arguments{
  \item{data}{A \eqn{n \times d}{n x d} matrix of component-wise maxima.}
  \item{x}{A \eqn{m \times d}{m x d} design matrix where the dependence function is evaluated, see \bold{Details}.}
  \item{d}{A postive integer (greater than or equal to two) indicating the number of variables. The trivariate case \code{d = 3} is the default.}
  \item{est}{A string denoting the estimation method (see \bold{Details}).}
    \item{margin}{A string denoting the type marginal distributions (see \bold{Details}).}
  \item{k}{A postive integer denoting the order of the Bernstein polynomial. \code{k = 13} is set by default.}
  \item{nboot}{A postive integer indicating the number of bootstrap replicates.}
  \item{y}{A numeric vector (of size \code{m}) with an initial estimate of the Pickands function. If \code{NULL}, The initial estimation is performed by using the estimation method chosen in \code{est}.}
  \item{matrix}{Logical; \code{FALSE} by default. If \code{TRUE}, and the dimension \code{d} is three, the value of \eqn{A} is collected in a square matrix.}
  \item{print}{Logical; \code{FALSE} by default. If \code{TRUE} the number of the iteration is printed.}
}

\details{
  Standard bootstrap is performed, in particular  estimates of the Pickands dependence function is provided for each data 
  resampling. 
  
  Most of the settings are the same as in the function \code{\link{beed}}.
    
  An empirical transformation of the marginals is performed when \code{margin = 'emp'}. Otherwise it refers to marginal parametric GEV estimation.
}

\value{
  \item{A}{Estimate of the Pickands function.}
  \item{bootA}{A matrix with \code{nboot} columns that reports the estimates of the Pickands function for each data resampling.}
  \item{beta}{A matrix of estimated polynomial coefficients. Each column is an estimates corrisponding to a data resampling.}
}

\references{
  Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. (2014)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials.
  \emph{Submitted}.
}

\seealso{\code{\link{beed}}, \code{\link{beed.confband}}.}

\examples{
# w <- seq(0 , 1, length = 100)
# data <- rbvevd(50, dep = 0.4, model = 'log', mar1 = c(1,1,1))
#  
# boot <- beed.boot(data, cbind(w, 1-w), 2, 'md', 'emp', 20, 500)  
#
}
\keyword{bootstrap}