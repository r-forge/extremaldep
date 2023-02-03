\name{returns}
\alias{returns}

\title{
  Compute return values
}

\description{
  Predicts the probability of future simultaneous exceedances
}

\usage{
returns(out, summary.mcmc, y, plot=FALSE, labels=NULL, 
        data=NULL)
}

\arguments{
  \item{out}{The output of the \code{\link{fExtDep.np}} function with \code{method="Bayesian"}.}
  \item{summary.mcmc}{The output of the \code{\link{summary_ExtDep}} function.}
  \item{y}{A 2-column matrix of unobserved thresholds.}
  \item{plot}{A logical value; if \code{TRUE}, then the returns are plotted using \code{\link{plot_ExtDep.np}}.}
  \item{labels}{As in \code{\link{plot_ExtDep.np}}. Required if \code{plot=TRUE}.}
  \item{data}{As in \code{\link{plot_ExtDep.np}}. Required if \code{plot=TRUE}.}    
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
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
	Giulia Marcon, \email{giuliamarcongm@gmail.com}
}

\references{
  Marcon, G., Padoan, S. A. and Antoniano-Villalobos, I. (2016).
  Bayesian inference for the extremal dependence.
  \emph{Electronic Journal of Statistics}, \bold{10}, 3310-3337.
}

\examples{

#########################################################
### Example 1 - daily log-returns between the GBP/USD ###
###             and GBP/JPY exchange rates            ###
#########################################################

if(interactive()){

data(logReturns)

mm_gbp_usd <- ts(logReturns$USD, start=c(1991,3), end=c(2014,12), frequency=12)
mm_gbp_jpy <- ts(logReturns$JPY, start=c(1991,3), end=c(2014,12), frequency=12)

### Detect seasonality and trend in the time series of maxima:
seas_usd <- stl(mm_gbp_usd, s.window="period")
seas_jpy <- stl(mm_gbp_jpy, s.window="period")

### remove the seasonality and trend from the two series:
mm_gbp_usd_filt <- mm_gbp_usd - rowSums(seas_usd$time.series[,-3])
mm_gbp_jpy_filt <- mm_gbp_jpy - rowSums(seas_jpy$time.series[,-3])

### Estimation of margins and dependence

mm_gbp <- cbind(as.vector(mm_gbp_usd_filt), as.vector(mm_gbp_jpy_filt))

hyperparam <- list(mu.nbinom = 3.2, var.nbinom = 4.48)
gbp_mar <- fExtDep.np(method="Bayesian", data=mm_gbp, par10=rep(0.1, 3), 
                      par20=rep(0.1,3), sig10=0.0001, sig20=0.0001, k0=5,
                      hyperparam = hyperparam, nsim=5e+4)

gbp_mar_sum <- summary_ExtDep(mcmc=gbp_mar, burn=3e+4, plot=TRUE)

mm_gbp_range <- apply(mm_gbp,2,quantile,c(0.9,0.995))

y_gbp_usd <- seq(from=mm_gbp_range[1,1], to=mm_gbp_range[2,1], length=20)
y_gbp_jpy <- seq(from=mm_gbp_range[1,2], to=mm_gbp_range[2,2], length=20)
y <- as.matrix(expand.grid(y_gbp_usd, y_gbp_jpy, KEEP.OUT.ATTRS = FALSE))

ret_marg <- returns(out=gbp_mar, summary.mcmc=gbp_mar_sum, y=y, plot=TRUE, 
                    data=mm_gbp, labels=c("GBP/USD exchange rate", "GBP/JPY exchange rate"))

}

#########################################################
### Example 2 - Reproducing some of the results shown ###
###             in Marcon et al. (2016, Figure 1)     ###
#########################################################

if(interactive()){

set.seed(1890)
data <- evd::rbvevd(n=100, dep=0.6, asy=c(0.8,0.3), model="alog", mar1=c(1,1,1))

hyperparam <- list(a.unif=0, b.unif=.5, mu.nbinom=3.2, var.nbinom=4.48)
pm0 <- list(p0=0.06573614, p1=0.3752118)

mcmc <- fExtDep.np(method="Bayesian", data=data, mar.fit=FALSE, k0=5,
                   pm0=pm0, prior.k = "nbinom", prior.pm = "unif", 
                   hyperparam=hyperparam, nsim=5e+5)

w <- seq(0.001, 0.999, length=100)
summary.mcmc <- summary_ExtDep(w, mcmc, burn=4e+5, plot=TRUE)

plot_ExtDep.np(out=mcmc, type = "A", summary.mcmc=summary.mcmc)
plot_ExtDep.np(out=mcmc, type = "h", summary.mcmc=summary.mcmc)
plot_ExtDep.np(out=mcmc, type = "pm", summary.mcmc=summary.mcmc)
plot_ExtDep.np(out=mcmc, type = "k", summary.mcmc=summary.mcmc)

y <- seq(10,100,2)
y <- as.matrix(expand.grid(y,y))
probs <- returns(out=mcmc, summary.mcmc=summary.mcmc, y=y, plot=TRUE)

}
}

\keyword{models}% use one of  RShowDoc("KEYWORDS")
