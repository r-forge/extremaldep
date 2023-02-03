\name{diagnostics}
\alias{diagnostics}

\title{Diagnostics plots for MCMC algorithm.}
\description{This function displays traceplots of the scaling parameter from the proposal distribution of the adaptive MCMC scheme and the associated acceptance probability.}

\usage{
diagnostics(mcmc)    
}
\arguments{
  \item{mcmc}{An output of the \code{\link{fGEV}} or \code{\link{fExtDep.np}} function with \code{method="Bayesian"}.}
}

\details{
  When \code{mcmc} is the output of \code{\link{fGEV}} then this corresponds to a marginal estimation and therefore \code{diagnostics} will  display in a first plot the value of \eqn{\tau}{\tau} the scaling parameter in the multivariate normal proposal which directly affects the acceptance rate of the proposal parameter values that are displayed in the second plot.
  
  When \code{mcmc} is the output of \code{\link{fExtDep.np}}, then this corresponds to an estimation of the dependence structure following the procedure given in Algorithm 1 of Beranger et al. (2021). If the margins are jointly estimated with the dependence (step 1 and 2 of the algorithm) then \code{diagnostics} provides trace plots of the corresponding scaling parameters (\eqn{\tau_1,\tau_2}{\tau_1,\tau_2}) and acceptance probabilities. For the dependence structure (step 3 of the algorithm), a trace plot of the polynomial order \eqn{\kappa}{\kappa} is given with the associated acceptance probability.
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
	Giulia Marcon, \email{giuliamarcongm@gmail.com}
}

\references{
	Beranger, B., Padoan, S. A. and Sisson, S. A. (2021).
	Estimation and uncertainty quantification for extreme quantile regions.
  \emph{Extremes}, \bold{24}, 349-375.
}

\seealso{\code{\link{fExtDep.np}}.}
\examples{

##################################################
### Example - Pollution levels in Milan, Italy ###
##################################################
	
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

f1 <- fGEV(data=data[,1], method="Bayesian", sig0 = 0.1, nsim = 5e+4)
diagnostics(f1)
burn1 <- 1:30000
gev.pars1 <- apply(f1$param_post[-burn1,],2,mean)
sdata1 <- trans2UFrechet(data=data[,1], pars=gev.pars1, type="GEV")

f2 <- fGEV(data=data[,2], method="Bayesian", sig0 = 0.1, nsim = 5e+4)
diagnostics(f2)
burn2 <- 1:30000
gev.pars2 <- apply(f2$param_post[-burn2,],2,mean)
sdata2 <- trans2UFrechet(data=data[,2], pars=gev.pars2, type="GEV")

sdata <- cbind(sdata1,sdata2)

### Bayesian estimation using Bernstein polynomials

pollut1 <- fExtDep.np(method="Bayesian", data=sdata, u=TRUE,
                      mar.fit=FALSE, k0=5, hyperparam = hyperparam, nsim=5e+4)

diagnostics(pollut1)

}
	
}

\keyword{Nonparametric}
