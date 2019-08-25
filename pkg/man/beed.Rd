\name{beed}
\alias{beed}
\title{Bernstein Polynomials Based Estimation of Extremal Dependence}
\description{
Estimates the Pickands dependence function corresponding to multivariate data on the basis of a Bernstein polynomials approximation.
}
\usage{
beed(data, x, d = 3, est = c("ht","cfg","md","pick"),
     margin = c("emp","est","exp","frechet","gumbel"),
     k = 13, y = NULL, beta = NULL, plot = FALSE)
}
\arguments{
  \item{data}{\eqn{(n \times d)}{(n x d)} matrix of component-wise maxima. \code{d} is the number of variables and \code{n} is the number of replications.}
  \item{x}{\eqn{(m \times d)}{(m x d)} design matrix where the dependence
    function is evaluated (see \bold{Details}).}
  \item{d}{positive integer greater than or equal to two indicating the number of variables (\code{d = 3} by default).}
  \item{est}{string, indicating the estimation method (\code{est="md"} by default, see \bold{Details}).}
  \item{margin}{string, denoting the type marginal distributions (\code{margin="emp"} by default, see \bold{Details}).}
  \item{k}{postive integer, indicating the order of Bernstein
    polynomials (\code{k=13} by default).}
  \item{y}{numeric vector (of size \code{m}) with an initial estimate of the Pickands function. If \code{NULL}, the initial estimate is computed using the estimation method denoted in \code{est}.}
  \item{beta}{vector of polynomial coefficients (see \bold{Details}).}
  \item{plot}{logical; if \code{TRUE} and \code{d<=3}, the estimated Pickands dependence function is plotted (\code{FALSE} by default).}
}


\details{
  The routine returns an estimate of the Pickands dependence function using the Bernstein polynomials approximation
  proposed in Marcon et al. (2017).
  The method is based on a preliminary empirical estimate of the Pickands dependence function.
  If you do not provide such an estimate, this is computed by the routine. In this case, you can select one of the empirical methods
  available. \code{est = 'ht'} refers to the Hall-Tajvidi estimator (Hall and Tajvidi 2000).
  With \code{est = 'cfg'} the method proposed by Caperaa et al. (1997) is considered. Note that in the multivariate case the adjusted  version of Gudendorf and Segers (2011) is used. Finally, with \code{est = 'md'} the estimate is based on the madogram defined in Marcon et al. (2017).

  Each row of the \eqn{(m \times d)}{(m x d)} design matrix \code{x} is a point in the unit \code{d}-dimensional simplex,

\eqn{
S_d := \left\{ (w_1,\ldots, w_d) \in [0,1]^{d}: \sum_{i=1}^{d} w_i = 1 \right\}.
}{
S_d := { (w_1,..., w_d) in [0,1]^{d}: \sum_{i=1}^{d} w_i = 1}.
}

  With this "regularization"" method, the final estimate satisfies the neccessary conditions in order to be a Pickands dependence function.

  \eqn{A(\bold{w}) = \sum_{\bold{\alpha} \in \Gamma_k} \beta_{\bold{\alpha}} b_{\bold{\alpha}} (\bold{w};k).}{A(w) = \sum_{\alpha in \Gamma_k} \beta_\alpha b_\alpha (w;k).}

  The estimates are obtained by solving an optimization quadratic problem subject to the constraints. The latter are represented
  by the following conditions:
  \eqn{A(e_i)=1; \max(w_i)\leq A(w) \leq 1; \forall i=1,\ldots,d;}{A(e_i)=1; max(w_i)<=A(w)<=1; for all i=1,...,d;} (convexity).

  The order of polynomial \code{k} controls the smoothness of the estimate. The higher \code{k} is, the smoother the final estimate is.
  Higher values are better with strong dependence (e. g. \code{k = 23}), whereas small values (e.g. \code{k = 6} or \code{k = 10}) are enough with mild or weak dependence.

An empirical transformation of the marginals is performed when \code{margin="emp"}. A max-likelihood fitting of the GEV distributions is implemented when \code{margin="est"}. Otherwise it refers to marginal parametric GEV theorethical distributions (\code{margin = "exp", "frechet", "gumbel"}).
}

\note{
The number of coefficients depends on both the order of polynomial \code{k} and the dimension \code{d}. The number of parameters is explained in Marcon et al. (2017).

The size of the vector \code{beta} must be compatible with the polynomial order \code{k} chosen.

With the estimated polynomial coefficients, the extremal coefficient, i.e. \eqn{d*A(1/d,\ldots,1/d)}{d*A(1/d,...,1/d)} is computed.
}

\value{
  \item{beta}{vector of polynomial coefficients}
  \item{A}{numeric vector of the estimated Pickands dependence function \eqn{A}}
  \item{Anonconvex}{preliminary non-convex function}
  \item{extind}{extremal index}
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.unibocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{http://www.borisberanger.com};
	Giulia Marcon, \email{giuliamarcongm@gmail.com}
}


\references{
  Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. and Segers, J. (2017)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials.
  \emph{Journal of Statistical Planning and Inference}, \bold{183}, 1-17.
}

\seealso{\code{\link{beed.confband}}.}

\examples{

x <- ExtremalDep:::simplex(2)
data <- evd::rbvevd(50, dep = 0.4, model = "log", mar1 = c(1,1,1))

Amd <- beed(data, x, 2, "md", "emp", 20, plot=TRUE)
Acfg <- beed(data, x, 2, "cfg", "emp", 20)
Aht <- beed(data, x, 2, "ht", "emp", 20)

lines(x[,1], Aht$A, lty = 1, col = 3)
lines(x[,1], Acfg$A, lty = 1, col = 2)

##################################
# Trivariate case
##################################

\donttest{

x <- ExtremalDep:::simplex(3)

data <- evd::rmvevd(50, dep = 0.8, model = "log", d = 3, mar = c(1,1,1))

par(mfrow=c(1,3))
Amd <- beed(data, x, 3, "md", "emp", 18, plot=TRUE)
Acfg <- beed(data, x, 3, "cfg", "emp", 18, plot=TRUE)
Aht <- beed(data, x, 3, "ht", "emp", 18, plot=TRUE)
par(mfrow=c(1,1))
}

}

\keyword{Nonparametric}
