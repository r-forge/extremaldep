\name{beed.boot}
\alias{beed.boot}
\title{Bootstrap Resampling and Bernstein Estimation of Extremal Dependence}
\description{
  Generate \code{nboot} bootstrap replicates from \code{data} and provide estimation of the uni- and multi-variate Pickands dependence function \eqn{A}, based on an Bernstein polynomials approximation, for each resample.
}
\usage{
  beed.boot(data, x, d = 3, est = c('ht','cfg','md'), margin = c(NULL, "exp-frechet", "frechet-exp"), k = 13, nboot = 500, y = NULL, beta = NULL, matrix = FALSE, print = FALSE)
}
\arguments{
  \item{data}{A matrix with \code{d} columns.}
  \item{x}{A \eqn{m \times d}{m x d} design matrix with \code{d} columns, in which the dependence function is evaluated, see \bold{Details}.}
  \item{d}{The dimension; an integer greater than or equal to two. The trivariate case \code{d = 3} is the default.}
  \item{est}{The estimation method (see \bold{Details}). Typically either \code{"cfg"} or \code{"ht"}. The method \code{"md"} is the default.}
    \item{margin}{Marginal transformation.}
  \item{k}{An integer parameter for the order of the Bernstein polynomials. \code{k = 13} is set by default.}
  \item{nboot}{An integer parameter for the number of bootstrap replicates.}
  \item{y}{Preliminary estimation. If \code{NULL}, The initial estimation is performed by using the estimation method chosen in \code{est}.}
  \item{beta}{A vector of beta coefficients (see \bold{Details}).}
  \item{matrix}{Logical; \code{FALSE} by default. If \code{TRUE}, and the dimension \code{d} is three, the value of \eqn{A} is collected in a square matrix.}
  \item{print}{Logical; \code{FALSE} by default. If \code{TRUE} the number of the iteration is printed.}
}

\details{
  Usual bootstrap is performed. For each resample, the estimation of the Pickands dependence function is provided. The setting is the same as in the function \code{\link{beed}}.
  
  For the \eqn{m \times d}{m x d} design matrix \code{x}, each row is a point in the unit \code{d}-dimensional simplex,

\eqn{
S_d := \left\{ (w_1,\ldots, w_{d-1}) \in [0,1]^{d-1}: \sum_{i=1}^{d-1} w_i \leq 1 \right\}.
}{
S_d := { (w_1,..., w_{d-1}) in [0,1]^{d-1}: \sum_{i=1}^{d-1} w_i <= 1}.
}
  
  The length of the vector \code{beta} must be compatible with the polynomial \code{k} chosen. The number of coefficients depends on both the order of polynomial \code{k} chosen and the dimension \code{d}. The relation is explicated in equation (8) in Marcon et al. (2014).
  
}

%\note{

\value{
  \item{A}{Estimate of the Pickands function of \code{data}.}
  \item{bootA}{A matrix with \code{nboot} columns which contains the estimates of the Pickands function of each resample.}
  \item{beta}{Matrix of coefficients \eqn{\beta}. Each column correspond to each resample.}
}

\references{
  Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. (2014)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials.
  \emph{Submitted}.
}

\seealso{\code{\link{beed}}, \code{\link{beed.confband}}.}

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
\keyword{bootstrap}






