\name{returns}
\alias{returns}

\title{
  Compute return values
}

\description{
  Predicts the probability of future simultaneous exceedances
}

\usage{
  returns(mcmc, summary.mcmc, y, plot=FALSE)
}

\arguments{
  \item{mcmc}{The output of the \code{bbeed} function.}
  \item{summary.mcmc}{The output of the \code{summary.bbeed} function.}
  \item{y}{A 2-column matrix of unobserved thresholds.}
  \item{plot}{If \code{plot=TRUE}, then the \code{plot.bbeed} function is used.}
}

\details{
  Computes for a range of unobserved extremes (larger than those observed in a sample), the pointwise mean from the posterior predictive distribution of such predictive values. The probabilities are calculated through

  \eqn{
P(Y_1 > y_1, Y_2 > y_2) = \frac{2}{k} \sum_{j=0}^{k-2} (\eta_{j+1} - \eta_j)
\times \left(
\frac{(j+1) B(y_1/(y_1+y_2)| j+2, k-j-1)}{y_1}
- \frac{(k-j-1) B(y_2/(y_1+y_2)| k-j, j+1)}{y_2}
\right),
}{P(Y_1 > y_1, Y_2 > y_2) = \frac{1}{k} \sum_{j=0}^{k-2} (\eta_{j+1} - \eta_j)
\times \left(
\frac{(j+1) B(y_1/(y_1+y_2)| j+2, k-j-1)}{y_1}
- \frac{(k-j-1) B(y_2/(y_1+y_2)| k-j, j+1)}{y_2}
\right),}

where \eqn{B(x|a,b)}{B(x|a,b)} denotes the cumulative distribution function of a Beta random variable with shape \eqn{a,b>0}{a,b>0}. See Marcon et al. (2016, p.3323) for details.
}

\value{
  Returns a vector whose length is equal to the number of rows of the input value \code{y}.
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

\examples{
		\dontrun{

	# This reproduces some of the results shown in Fig. 1 (Marcon, 2016).
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

	y <- seq(10,100,2)
	y <- as.matrix(expand.grid(y,y))
	probs <- returns(mcmc = mcmc, summary.mcmc = summary.mcmc, y = y, plot = TRUE)

	}
}

\keyword{ models }% use one of  RShowDoc("KEYWORDS")
