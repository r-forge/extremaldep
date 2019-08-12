\name{plot.bbeed}
\alias{plot.bbeed}
\title{Plot of Extremal Dependence}
\description{
  Produces one or more plots of the extremal dependence.
}
\usage{
\method{plot}{bbeed}(x, type = c("summary", "returns", "A", "h", "pm", "k"),
          mcmc, summary.mcmc, nsim, burn, y, probs, CEX=1.5, A_true, h_true,
          labels=c(expression(y[1]),expression(y[2])), ...)
}
\arguments{
  \item{x}{Vector on the unit simplex where the dependence function is evaluated.}
  \item{type}{String, denoting the type of function to plot (see \bold{Details}).}
  \item{mcmc}{The output of the bbeed function.}
  \item{summary.mcmc}{The output of the summary.bbeed function.}
  \item{nsim}{The number of simulation in the mcmc algorithm.}
  \item{burn}{The burn-in period.}
  \item{y}{A 2-column matrix of unobserved thresholds at which the returns are calculated. Required when \code{type="y"}.}
  \item{probs}{The probability of joint exceedances, the output of the \code{return} function.}
  \item{A_true}{The true pickands dependence function (evaluated at \code{x}).}
  \item{h_true}{The true angular density function (evaluated at \code{x}).}
  \item{CEX}{Label and axis sizes.}
  \item{labels}{Labels.}
  \item{...}{Additional graphical parameters. See \code{plot} function for details.}
}

\details{
	If \code{type="returns"}, a (contour) plot of the probabilities of exceedances for some threshold \code{y} is returned. \cr
	If \code{type="A"}, a plot of the estimated Pickands dependence function is drawn. If \code{A_true} is specified the plot includes the true Pickands dependence function and a functional boxplot for the estimated function.
	If \code{type="h"}, a plot of the estimated angular density function is drawn. If \code{h_true} is specified the plot includes the true angular density and a functional boxplot for the estimated function.
	If \code{type="pm"}, a plot of the prior against the posterior for the mass at \eqn{\{0\}}{\{0\}} is drawn.
	If \code{type="k"}, a plot of the prior against the posterior for the polynomial degree \eqn{k}{k} is drawn.
	If \code{type="summary"}, a 2 by 2 plot with types \code{"A"}, \code{"h"}, \code{"pm"} and \code{"k"} is returned.
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.unibocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{http://www.borisberanger.com};
	Giulia Marcon, \email{giuliamarcongm@gmail.com}
}

\references{
  Marcon, G., Padoan, S.A., Naveau, P., Muliere, P., Segers, J. (2017)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials.
  \emph{Journal of Statistical Planning and Inference}, \bold{183}, 1-17.
}

\seealso{\code{\link{simplex}}, \code{\link{beed.confband}}.}
\examples{
	\dontrun{

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

	plot.bbeed(type = "A", x=w, mcmc=mcmc, summary.mcmc, nsim=nsim, burn=burn)
	plot.bbeed(type = "h", x=w, mcmc=mcmc, summary.mcmc, nsim=nsim, burn=burn)
	plot.bbeed(type = "pm", x=w, mcmc=mcmc, summary.mcmc, nsim=nsim, burn=burn)
	plot.bbeed(type = "k", x=w, mcmc=mcmc, summary.mcmc, nsim=nsim, burn=burn)

	Atrue <- evd::abvevd(w, dep=0.6, asy=c(0.3,0.8), model='alog')
	htrue <- evd::hbvevd(w, dep=0.6, asy=c(0.8,0.3), model='alog',half=TRUE)

	par(mfrow=c(1,1))
	plot.bbeed(type = "A", summary.mcmc=summary.mcmc, A_true=Atrue)
	plot.bbeed(type = "h", summary.mcmc=summary.mcmc, h_true=htrue)

	}
}

\keyword{Nonparametric}
