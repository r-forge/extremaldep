\name{beed}
\alias{beed}
\title{Bernstein Estimation of Extremal Dependence}
\description{
  Estimate an uni- or multi-variate Pickands dependence function
  on the basis of a Bernstein polynomials approximation
}
\usage{
beed(data, x, d = 3, est = c('ht','cfg','md'), margin = c('emp','Gev'),
     k = 13, y = NULL, beta = NULL, matrix = FALSE, plot = FALSE)
}
\arguments{
  \item{data}{A \eqn{(n \times d)}{(n x d)} matrix of component-wise maxima. \code{d} is the number of variables and \code{n} is the number of replications.}
  \item{x}{A \eqn{(m \times d)}{(m x d)} design matrix where the dependence
    function is evaluated (see \bold{Details}).}
  \item{d}{A positive integer greater than or equal to two indicating the number of variables (\code{d = 3} is the default).}
  \item{est}{A string, indicating the estimation method (see \bold{Details}). 
  The method \code{'md'} is the default.}
  \item{margin}{A string denoting the type marginal distributions (see \bold{Details}).}
  \item{k}{A postive integer indicating the order of the Bernstein 
    polynomials. The value \code{k = 13} is the default.}
  \item{y}{A numeric vector (of size \code{m}) with an initial estimate of the Pickands function. If \code{NULL}, the initial 
    estimate is computed using the estimation method denoted in
    \code{est}.}
  \item{beta}{A vector of polynomial coefficients (see \bold{Details}).}
  \item{matrix}{Logical; \code{FALSE} by default. If \code{TRUE}, and the dimension \code{d} is three (the default dimension), the value of
    \eqn{A} is collected in a square matrix.}
  \item{plot}{Logical; if \code{TRUE} and \code{d=3} or \code{d=2}
    the estimated Pickands dependence function is plotted. If \code{FALSE}
    (the default), the plot is not provided.}
}

\details{
  The routine returns an estimate of the Pickands dependence function using the Bernstein polynomials approximation 
  proposed in Marcon et al. (2014). 
  The method is based on a preliminary empirical estimate of the Pickands dependence function. 
  If you do not provide such an estimate, this is computed by the routine. In this case, you can select one of empirical methods
  available. \code{est = 'ht'} refers to the Hall-Tajvidi estimator (Hall and Tajvidi 2000). 
  With \code{est = 'cfg'} the method proosed by Caperaa et al. (1997) is considered. Note that in the multivariate case the adjusted  version of Gudendorf and Segers (2011) is used. Finally, with \code{est = 'md'} the estimate is based on the madogram defined in Marcon et al. (2014). 
  
  Each row of the \eqn{(m \times d)}{(m x d)} design matrix \code{x} is a point in the unit \code{d}-dimensional simplex,
  
\eqn{
S_d := \left\{ (w_1,\ldots, w_{d-1}) \in [0,1]^{d-1}: \sum_{i=1}^{d-1} w_i \leq 1 \right\}.
}{
S_d := { (w_1,..., w_{d-1}) in [0,1]^{d-1}: \sum_{i=1}^{d-1} w_i <= 1}.
}

  With this 'regularization' method, the final estimate satisfies the neccessary conditions in order to be a Pickands dependence function. 
  
  \eqn{A(\bold{w}) = sum_j \beta_j b(j,\bold{w};k)}{A(w) = sum_j beta_j b(j,w;k)}
  
  The estimates are obtained solving an optimization quadratic problem subject to the constraints. The latter are represented 
  by the following conditions:
  \eqn{A(e_i)=1; max(w_i)\leq A(w) \leq 1; \forall i=1,\ldots,d; convexity.}{A(e_i)=1; max(w_i)<=A(w)<=1; for all i=1,...,d; convexity.}
  
  The order of polynomial \code{k} controls the smoothness of the estimate. The higher \code{k} is, the smoother the final estimate is.
  Higher values are better with strong dependence (e. g. \code{k = 23}), whereas small values (e.g. \code{k = 6} or \code{k = 10}) are enough with mild or weak dependence.
  
  An empirical transformation of the marginals is performed when \code{margin = 'emp'}. Otherwise it refers to marginal parametric GEV estimation.
}

\note{
The number of coefficients depends on both the order of polynomial \code{k} and the dimension \code{d}. The number of parameters is explained in Marcon et al. (2014).

The size of the vector \code{beta} must be compatible with the polynomial order \code{k} chosen.

With the estimated polynomial coefficients, the extremal coefficient, i.e. \eqn{d*A(1/d,\ldots,1/d)}{d*A(1/d,...,1/d)} is computed. 
}

\value{
  \item{beta}{vector of polynomial coefficients}
  \item{A}{the Pickands dependence function \eqn{A} estimated}
  \item{Anonconvex}{preliminary non-convex function}
  \item{extind}{The extremal index}
}

\author{
  Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{faculty.bocconi.it/simonepadoan};
	Giulia Marcon
}

\references{
  Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. and Segers, J. (2016)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials.
  \emph{Journal of Statistical Planning and Inference}, \bold{To appear}.
}

\seealso{\code{\link{beed.confband}}.}

\examples{

w <- seq(0, 1, length = 100)  
data <- rbvevd(50, dep = 0.4, model = 'log', mar1 = c(1,1,1))
    
Amd <- beed(data, cbind(w, 1-w), 2, 'md', 'emp', 20, plot=TRUE) 
Acfg <- beed(data, cbind(w, 1-w), 2, 'cfg', 'emp', 20) 
Aht <- beed(data, cbind(w, 1-w), 2, 'ht', 'emp', 20) 

lines(w, Aht$A, lty = 1, col = 3)
lines(w, Acfg$A, lty = 1, col = 2)

##################################
# Trivariate case
##################################

\dontrun{
w <- seq(0, 1, length = 30)
xx <- expand.grid(w, w)
xx <- cbind(xx, 1 - rowSums(xx))
os <- rowSums(xx[,1:2])>1
xx[os,] <- rep(NA,3)
x <- as.matrix(xx)

data <- rmvevd(50, dep = 0.8, model = "log", d = 3, mar = c(1,1,1))

par(mfrow=c(1,3))
Amd <- beed(data, x, 3, 'md', 'emp', 18, plot=TRUE)
Acfg <- beed(data, x, 3, 'cfg', 'emp', 18, plot=TRUE)
Aht <- beed(data, x, 3, 'ht', 'emp', 18, plot=TRUE)
}

}

\keyword{Nonparametric}