\name{beed}
\alias{beed}
\title{Bernstein Estimation of Extremal Dependence}
\description{
  Estimate of the uni- and multi-variate Pickands dependence function
  \eqn{A} based on an Bernstein polynomials approximation
}
\usage{
beed(data, x, d = 3, est = c('ht','cfg','md'), k = 13, y = NULL, 
    beta = NULL, matrix = FALSE, plot = FALSE)
}
\arguments{
  \item{data}{A matrix with \code{d} columns, which
    contain point wise maxima.}
  \item{x}{A \eqn{m \times d}{m x d} design matrix (see \bold{Details}) in which the dependence
    function is evaluated.}
  \item{d}{The dimension; an integer greater than or equal to two.
    The trivariate case \code{d = 3} is the default.}
  \item{est}{The estimation method (see \bold{Details}). Typically
    either \code{"cfg"} or \code{"ht"}. The method \code{"md"} is the
    default.}
  \item{margin}{Marginal transformation.}
  \item{k}{An integer parameter for the order of the Bernstein 
    polynomials. \code{k = 13} is set by default.}
  \item{y}{Preliminary estimation. If \code{NULL}, The initial 
    estimation is performed by using the estimation method chosen in
    \code{est}.}
  \item{beta}{A vector of beta coefficients (see \bold{Details}). The
    length must be compatible with the polynomial order \code{k} chosen.}
  \item{matrix}{Logical; \code{FALSE} by default. If \code{TRUE}, and the dimension \code{d} is three (the default dimension), the value of
    \eqn{A} is collected in a square matrix.}
  \item{plot}{Logical; if \code{TRUE}, and the dimension \code{d} is
    three (the default dimension), the contours of the dependence function of a trivariate extreme value distribution are plotted. For plotting in the bivariate case, use \code{d=2}. If \code{FALSE}
    (the default), the following arguments are ignored.}
}

\details{
  The procedure recalls the estimation methodology proposed in Marcon et al. (2014) and it is based on the Bernstein polynomials approximation. 
  A preliminary estimation is required. If it's not available, it's computed specifying the method preferred among those proposed in literature. \code{est = "ht"} refers to Hall and Tajvidi (2000); \code{est = "cfg"} use the method of Caperaa et al. (1997) for the bivariate case and the adjusted one of Gudendorf and Segers (2011) for the multivariate case; \code{est = "md"} is based on the madogram defined in Marcon et al. (2014). 
  For the \eqn{m \times d}{m x d} design matrix \code{x}, each row is a point in the unit \code{d}-dimensional simplex,

\eqn{
S_d := \left\{ (w_1,\ldots, w_{d-1}) \in [0,1]^{d-1}: \sum_{i=1}^{d-1} w_i \leq 1 \right\}.
}{
S_d := { (w_1,..., w_{d-1}) in [0,1]^{d-1}: \sum_{i=1}^{d-1} w_i <= 1}.
}
  The final estimation is the Pickands dependence function \eqn{A} where all necessary constraints are satisfied. 
  
  \eqn{A(\bold{w}) = sum_j \beta_j b(j,\bold{w};k)}{A(w) = sum_j beta_j b(j,w;k)}
  
  The routine is based on an optimization quadratic problem subject to the constraints required and identified by some particular matrices:
  \eqn{A(e_i)=1; max(w_i)\leq A(w) \leq 1; \forall i=1,\ldots,d; convexity.}{A(e_i)=1; max(w_i)<=A(w)<=1; for all i=1,...,d; convexity.}
  
  Setting the order of polynomial \code{k} allows to chose how smooth the final estimation should be. Higher values are preferable when the dependence is strong (e. g. \code{k = 23}), while values such as \code{k = 10} or \code{k = 6} are enough with mild and weak dependence, respectively.
}

\note{
  Once the vector of coefficients \code{beta} is estimated, it can be used to compute the extremal coefficient as \eqn{d*A(1/d,\ldots,1/d)}{d*A(1/d,...,1/d)}. Notice that the number of coefficients depends on both the order of polynomial \code{k} chosen and the dimension \code{d}. The relation is explicated in equation (8) in Marcon et al. (2014).
}

\value{
  A list with two elements corresponding to the Pickands dependence
  function \eqn{A} estimated and the vector of coefficient \eqn{\beta}{beta}.
}

\references{
  Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. (2014)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials.
  \emph{Submitted}.
}

\seealso{\code{\link{beed.confband}}.}

\examples{
w <- seq(0,1, length=100)  
data <- rbvevd(50, dep = 0.4, model = "log", mar1 = c(1,1,1))
    
Amd <- beed(data = data, x = cbind(w, 1-w), d = 2, k = 20, 
            est = 'md', margin = "frechet-exp")
Acfg <- beed(data = data, x = cbind(w, 1-w), d = 2, k = 20, 
             est = 'md', margin = NULL)
Aht <- beed(data = data, x = cbind(w, 1-w), d = 2, k = 20, 
            est = 'md', margin = NULL)

plot(w, w, type = "n", xlab = "w1", ylab = "w2", ylim = c(.5,1))
polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = 1, lwd = 1, border = 'grey')
lines(w, Aht$A, lty = 1, col = 3)
lines(w, Acfg$A, lty=1, col = 2)
lines(w, Amd$A, lty=1, col = 1)
}

\keyword{nonparametric}







