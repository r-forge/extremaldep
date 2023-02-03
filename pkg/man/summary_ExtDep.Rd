\name{summary_ExtDep}
\alias{summary_ExtDep}

\title{Summary of MCMC algorithm.}
\description{This function computes summaries on the posterior sample obtained from the adaptive MCMC scheme for the non-parametric estimation of a bivariate dependence structure.}

\usage{
summary_ExtDep(object, mcmc, burn, cred=0.95, plot=FALSE, ...)    
}
\arguments{
  \item{object}{A vector of values on \eqn{[0,1]}{[0,1]}. If missing, a regular grid of length \eqn{100}{100} is considered.}
  \item{mcmc}{An output of the \code{\link{fExtDep.np}} function with \code{method="Bayesian"}.}
  \item{burn}{A positive integer indicating the burn-in period.}
  \item{cred}{A value in \eqn{[0,1]}{[0,1]} indicating the level of the credibility intervals to be computed.}
  \item{plot}{A logical value; if \code{TRUE} a summary of the estimated dependence is displayed by calling the \code{plot_ExtDep.np} function with \code{type="summary"}.}
  \item{...}{Additional graphical parameters for \code{plot_ExtDep.np} when \code{plot=TRUE}.}
}

\details{

  For each value say \eqn{\omega \in [0,1]}{\omega \in [0,1]} given, the complement \eqn{1-\omega }{1-\omega} is automatically computed to define the observation \eqn{(\omega,1-\omega)}{(\omega,1-\omega)} on the bivariate unit simplex.

  It is obvious that the value of \code{burn} must be greater than the number of iterations in the mcmc algorithm. This can be found in \code{mcmc}.
}

\value{
  The function returns a list with the following objects:
  \itemize{
  \item{\code{k.median}, \code{k.up}, \code{k.low}: }{Posterior median, upper and lower bounds of the CI for the estimated Bernstein polynomial degree \eqn{\kappa}{\kappa};}
  \item{\code{h.mean}, \code{h.up}, \code{h.low}: }{Posterior mean, upper and lower bounds of the CI for the estimated angular density \eqn{h}{h};}
  \item{\code{A.mean}, \code{A.up}, \code{A.low}: }{Posterior mean, upper and lower bounds of the CI for the estimated Pickands dependence function \eqn{A}{A};}
  \item{\code{p0.mean}, \code{p0.up}, \code{p0.low}: }{Posterior mean, upper and lower bounds of the CI for the estimated point mass \eqn{p_0}{p_0};}
  \item{\code{p1.mean}, \code{p1.up}, \code{p1.low}: }{Posterior mean, upper and lower bounds of the CI for the estimated point mass \eqn{p_1}{p_1};}
  \item{\code{A_post}: }{Posterior sample for Pickands dependence function;}
  \item{\code{h_post}: }{Posterior sample for angular density;}
  \item{\code{eta.diff_post}: }{Posterior sample for the Bernstein polynomial coefficients (\eqn{\eta}{\eta} parametrisation);}
  \item{\code{beta_post}: }{Posterior sample for the Bernstein polynomial coefficients (\eqn{\beta}{\beta} parametrisation);}
  \item{\code{p0_post}, \code{p1_post}: }{Posterior sample for point masses \eqn{p_0}{p_0} and \eqn{p_1}{p_1};}
  \item{\code{w}: }{A vector of values on the bivariate simplex where the angular density and Pickands dependence function were evaluated;}
  \item{\code{burn}: }{The argument provided;}
  }
  If the margins were also fitted, the list given as \code{object} would contain \code{mar1} and \code{mar2} and the function would also output:
  \itemize{
  \item{\code{mar1.mean}, \code{mar1.up}, \code{mar1.low}: }{Posterior mean, upper and lower bounds of the CI for the estimated marginal parameter on the first component;}
  \item{\code{mar2.mean}, \code{mar2.up}, \code{mar2.low}: }{Posterior mean, upper and lower bounds of the CI for the estimated marginal parameter on the second component;}
  \item{\code{mar1_post}: }{Posterior sample for the estimated marginal parameter on the first component;}
  \item{\code{mar2_post}: }{Posterior sample for the estimated marginal parameter on the second component;}  
  }
}


\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com}
}

\seealso{\code{\link{fExtDep.np}}.}
\examples{

####################################################
### Example - Pollution levels in Milan, Italy ###
####################################################
	
if(interactive()){

### Here we will only model the dependence structure	
data(MilanPollution)

data <- Milan.winter[,c("NO2","SO2")] 
data <- as.matrix(data[complete.cases(data),])

# Thereshold
u <- apply(data, 2, function(x) quantile(x, prob=0.9, type=3))

# Hyperparameters
hyperparam <- list(mu.nbinom = 6, var.nbinom = 8, a.unif=0, b.unif=0.2)

### Standardise data to univariate Frechet margins

f1 <- fGEV(data=data[,1], method="Bayesian", sig0 = 0.0001, nsim = 5e+4)
diagnostics(f1)
burn1 <- 1:30000
gev.pars1 <- apply(f1$param_post[-burn1,],2,mean)
sdata1 <- trans2UFrechet(data=data[,1], pars=gev.pars1, type="GEV")

f2 <- fGEV(data=data[,2], method="Bayesian", sig0 = 0.0001, nsim = 5e+4)
diagnostics(f2)
burn2 <- 1:30000
gev.pars2 <- apply(f2$param_post[-burn2,],2,mean)
sdata2 <- trans2UFrechet(data=data[,2], pars=gev.pars2, type="GEV")

sdata <- cbind(sdata1,sdata2)

### Bayesian estimation using Bernstein polynomials

pollut1 <- fExtDep.np(method="Bayesian", data=sdata, u=TRUE,
                      mar.fit=FALSE, k0=5, hyperparam = hyperparam, nsim=5e+4)

diagnostics(pollut1)
pollut1_sum <- summary_ExtDep(mcmc=pollut1, burn=3e+4, plot=TRUE)
}
	
}

\keyword{Nonparametric}
