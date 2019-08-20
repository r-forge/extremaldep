\name{beed.boot}
\alias{beed.boot}
\title{Bootstrap Resampling and Bernstein Estimation of Extremal Dependence}
\description{
 Computes \code{nboot} estimates of the Pickands dependence function for multivariate data (using the Bernstein polynomials approximation method) on the basis of the bootstrap resampling of the data.
}
\usage{
  beed.boot(data, x, d = 3, est = c("ht","md","cfg","pick"),
     margin=c("emp", "est", "exp", "frechet", "gumbel"), k = 13,
     nboot = 500, y = NULL, print = FALSE)
}
\arguments{
  \item{data}{\eqn{n \times d}{n x d} matrix of component-wise maxima.}
  \item{x}{\eqn{m \times d}{m x d} design matrix where the dependence function is evaluated, see \bold{Details}.}
  \item{d}{postive integer (greater than or equal to two) indicating the number of variables (\code{d=3} by default).}
  \item{est}{string denoting the preliminary estimation method (see \bold{Details}).}
    \item{margin}{string denoting the type marginal distributions (see \bold{Details}).}
  \item{k}{postive integer denoting the order of the Bernstein polynomial (\code{k=13} by default).}
  \item{nboot}{postive integer indicating the number of bootstrap replicates (\code{nboot=500} by default).}
  \item{y}{numeric vector (of size \code{m}) with an initial estimate of the Pickands function. If \code{NULL}, The initial estimation is performed by using the estimation method chosen in \code{est}.}
  \item{print}{logical; \code{FALSE} by default. If \code{TRUE} the number of the iteration is printed.}
}

\details{
  Standard bootstrap is performed, in particular estimates of the Pickands dependence function are provided for each data resampling.

  Most of the settings are the same as in the function \code{\link{beed}}.

  An empirical transformation of the marginals is performed when \code{margin="emp"}. A max-likelihood fitting of the GEV distributions is implemented when \code{margin="est"}. Otherwise it refers to marginal parametric GEV theorethical distributions (\code{margin = "exp", "frechet", "gumbel"}).
}

\value{
  \item{A}{numeric vector of the estimated Pickands dependence function.}
  \item{bootA}{matrix with \code{nboot} columns that reports the estimates of the Pickands function for each data resampling.}
  \item{beta}{matrix of estimated polynomial coefficients. Each column corresponds to a data resampling.}
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
  \emph{Journal of Statistical Planning and Inference}, 183, 1-17.
}

\seealso{\code{\link{beed}}, \code{\link{beed.confband}}.}

\examples{
\dontrun{

x <- ExtremalDep:::simplex(2)
data <- rbvevd(50, dep = 0.4, model = "log", mar1 = c(1,1,1))

boot <- beed.boot(data, x, 2, "md", "emp", 20, 500)
}
}

\keyword{bootstrap}
