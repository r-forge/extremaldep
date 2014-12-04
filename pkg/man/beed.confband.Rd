\name{beed.confband}
\alias{beed.confband}
\title{Nonparametric Bootstrap Confidence Intervals}
\description{
  This function generates 2 different types of \eqn{(1-\alpha)\%} nonparametric bootstrap confident bands for the estimated Pickands dependence function.
}
\usage{
  beed.confband(data, x, d = 3, est = c('ht','cfg','md'), margin = c(NULL, "exp-frechet", "frechet-exp"), k = 13, nboot = 500, conf = 0.95, y = NULL, beta = NULL, matrix = FALSE, print = FALSE)
}
\arguments{
  \item{data}{A matrix with \code{d} columns}
  \item{x}{A \eqn{m \times d}{m x d} design matrix (see \bold{Details}).}
  \item{d}{The dimension; an integer greater than or equal to two. The trivariate case \code{d = 3} is the default.}
  \item{est}{The estimation method (see \bold{Details}). Typically either \code{"cfg"} or \code{"ht"}. The method \code{"md"} is the default.}
    \item{margin}{Marginal transformation.}
  \item{k}{An integer parameter for the order of the Bernstein polynomials. \code{k = 13} is set by default.}
  \item{nboot}{An integer parameter for the number of bootstrap replicates.}
  \item{conf}{A scalar containing the confidence level of the required interval. \code{conf = 0.95}{is by default.}}
  \item{y}{Preliminary estimation. If \code{NULL}, The initial estimation is performed by using the estimation method chosen in \code{est}.}
  \item{beta}{A vector of beta coefficients (see \bold{Details}). The length must be compatible with the polynomial \code{k} chosen.}
  \item{matrix}{Logical; \code{FALSE} by default. If \code{TRUE}, and the dimension \code{d} is three (the default dimension), the value of \eqn{A} is collected in a square matrix.}
  \item{print}{Logical; \code{FALSE} by default. If \code{TRUE}, the number of the iteration is printed.}
}

\details{
For the \eqn{m \times d}{m x d} design matrix \code{x}, each row is a point in the unit \code{d}-dimensional simplex,

\eqn{
S_d := \left\{ (w_1,\ldots, w_{d-1}) \in [0,1]^{d-1}: \sum_{i=1}^{d-1} w_i \leq 1 \right\}.
}{
S_d := { (w_1,..., w_{d-1}) in [0,1]^{d-1}: \sum_{i=1}^{d-1} w_i <= 1}.
}
Compute \eqn{(1-\alpha)\%} bootstrap confident bands for the estimated Pickands dependence function based on Bernstein polynomial representation. It uses the usual bootstrap routine (see \code{\link{beed.boot}}).
Two methods are provided. The former takes into account the quantiles of the bootstrap sample distribution of the estimated Pickands functions and computes the confidence bands pointwise. The latter is based on the bootstrap sampling distribution of the \eqn{\beta} coefficients. Once the quantiles are detected, they are used as coefficients in the Bernstein representation to construct the curves that contain the true Pickands function with, at most, \eqn{(1-\alpha)} of probability.

  The length of the vector \code{beta} must be compatible with the polynomial \code{k} chosen. The number of coefficients depends on both the order of polynomial \code{k} chosen and the dimension \code{d}. The relation is explicated in equation (8) in Marcon et al. (2014).
}

%\note{
  
  \value{
    \item{A}{Estimate of the Pickands function of \code{data}.}
    \item{bootA}{A matrix with \code{nboot} columns which contains the estimates of the Pickands function of each resample.}
    \item{A.up.beta/A.low.beta}{Vectors of upper and lower bands of the Pickands dependence function obtained considering the bootstrap sampling distribution of the \eqn{\beta} coefficients.}
    \item{A.up.pointwise/A.low.pointwise}{Vectors of upper and lower bands of the Pickands dependence function obtained considering the bootstrap sampling distribution of the estimated Pickands functions.}
    \item{up.beta/low.beta}{Vectors of upper and lower bounds of the bootstrap sampling distribution of the \eqn{\beta} coefficients.}
  }
  
  \references{
    Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. (2014)
    Multivariate Nonparametric Estimation of the Pickands   
    Dependence Function using Bernstein Polynomials.
    \emph{Submitted}.
  }
  
  \seealso{\code{\link{beed}}, \code{\link{beed.boot}}.}
  

\examples{
w <- seq(0,1, length=100)  
data <- rbvevd(50, dep = 0.4, model = "log", mar1 = c(1,1,1))
  
cb <- beed.confband(data = data, x = cbind(w, 1-w), d = 2, 
                    k = 20, nboot = 500, est = 'md', 
                    margin = "frechet-exp")

plot(w, w, type = 'n', xlab = 'w1', ylab = 'w2', ylim = c(.5,1))
polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = 1, lwd = 1, border = 'grey')
lines(w,cb$A,lty=1)
lines(w,cb$A.up.beta,lty=2)
lines(w,cb$A.low.beta,lty=2)
}
  
  \keyword{bootstrap confidence bands}
  
  
  
  
  
  