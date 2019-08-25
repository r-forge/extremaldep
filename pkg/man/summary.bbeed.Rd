\name{summary.bbeed}
\alias{summary.bbeed}

\title{
  Compute summary statistics from the MCMC output.
}
\description{
  Summary statistics of the MCMC output obtained from the Bayesian method based on the Bernstein polynomials for inferring the angular measure and Pickands dependence functions.
}
\usage{
  \method{summary}{bbeed}(object, mcmc, burn, conf=0.95, plot=FALSE, ...)
}

\arguments{
  \item{object}{The values (on the simplex) at which the dependence is evaluated}
  \item{mcmc}{The output of an MCMC algorithm given by \code{beed}}
  \item{burn}{The burn-in period}
  \item{conf}{The confidence region}
  \item{plot}{If \code{plot=TRUE} then the function \code{plot.beed} is used.}
  \item{...}{Arguments to be passed for the graphical parameters}
}

\value{
	Returns a list that contains:
	\cr
	* the \code{conf}-, 0.5- and \code{1-conf}-quantiles and posterior sample for k (the polynomial order),
	\cr
	* the \code{conf}- and \code{1-conf}-quantiles, mean and posterior sample for h (the angular density), A (the Pickands dependence function), p0 and p1 (the point masses at the endpoints of the simplex), mar1 and mar2 (the marginal parameters, if they exist). To access them, the names are for example \code{k.low}, \code{k.median}, \code{k.up} and \code{k_post}.
	\cr
	* \code{w} and \code{burn} which are the inputs \code{object} and \code{burn}.
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.unibocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{http://www.borisberanger.com};
	Giulia Marcon, \email{giuliamarcongm@gmail.com}
}

\references{
  Marcon G., Padoan, S.A. and Antoniano-Villalobos I. (2016) Bayesian Inference for the Extremal Dependence. \emph{Electronic Journal of Statistics}, 10.2, 3310-3337.
}

\seealso{\code{\link{plot.bbeed}}.}

\examples{
\donttest{

# This reproduces some of the results showed in Fig. 1 (Marcon, 2016).
set.seed(1890)
data <- evd::rbvevd(n=100, dep=.6, asy=c(0.8,0.3), model="alog", mar1=c(1,1,1))

nsim = 500000
burn = 400000

mu.nbinom = 3.2
var.nbinom = 4.48
hyperparam <- list(a.unif=0, b.unif=.5, mu.nbinom=mu.nbinom, var.nbinom=var.nbinom)
k0 = 5
pm0 = list(p0=0.06573614, p1=0.3752118)
eta0 = ExtremalDep:::rcoef(k0, pm0)

mcmc <- bbeed(data, pm0, eta0, k0, hyperparam, nsim,
              prior.k = "nbinom", prior.pm = "unif")

w <- seq(0.001, .999, length=100)
summary.mcmc <- summary.bbeed(w, mcmc, burn, nsim, plot=TRUE)
}
}

\keyword{aplot}% use one of  RShowDoc("KEYWORDS")
