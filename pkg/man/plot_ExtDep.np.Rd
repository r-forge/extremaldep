\name{plot_ExtDep.np}
\alias{plot_ExtDep.np}

\title{Graphical summaries of non-parametric representations of extremal dependence.}
\description{This function displays several summaries of extremal dependence represented through Bernstein polynomials.}

\usage{
plot_ExtDep.np(out, type, summary.mcmc, burn, y, probs, 
               A_true, h_true, est.out, mar1, mar2, dep, 
               QatCov1=NULL, QatCov2=QatCov1, P, 
               labels=c(expression(y[1]),expression(y[2])), 
               CEX=1.5, xlim, ylim, col.data, 
               col.Qfull, col.Qfade, data=NULL, ...)        
}
\arguments{
  \item{out}{An output of the \code{\link{fExtDep.np}} function.}
  \item{type}{A character string indicating the type of graphical summary to be plotted. Takes values \code{"summary"}, \code{"returns"}, \code{"A"}, \code{"h"}, \code{"pm"}, \code{"k"} or \code{"Qsets"}.}
  \item{summary.mcmc}{The output of the \code{summary_ExtDep} function. Only required when \code{out} is obtained using a Bayesian estimation method (\code{out$est=="Bayesian"}).}
  \item{burn}{The burn-in period. Only required when \code{out} is obtained using a Bayesian estimation method (\code{out$est=="Bayesian"}).}
  \item{y}{A 2-column matrix of unobserved thresholds at which the returns are calculated. Required when \code{type="returns"}.}
  \item{probs}{The probability of joint exceedances, the output of the \code{\link{return}} function.}
  \item{A_true}{A vector representing the true pickands dependence function evaluated at the grid points on the simplex given by \code{summary.mcmc$w}.}
  \item{h_true}{A vector representing the true angular density function evaluated at the grid points on the simplex given by \code{summary.mcmc$w}.}
  \item{est.out}{A list containing:
    \itemize{
      \item{\code{ghat}: }{a 3-row matrix giving the posterior summary for the estimate of the bound;}
      \item{\code{Shat} and \code{Shat_post}: }{a matrix of the posterior and a 3-row matrix giving the posterior summary for the estimate of the basic set \eqn{S}{S};}
      \item{\code{nuShat} and \code{nuShat_post}: }{a matrix of the posterior and a 3-row matrix giving the posterior summary for the estimate of the measure \eqn{\nu(S)}{\nu(S)};}
    }
    Note that a posterior summary is made of its mean and \eqn{90\%}{90\%}credibility interval.
    
    Only required when using a Bayesian estimation method (\code{out$est=="Bayesian"}) and \code{type="Qsets"}.}
  \item{mar1, mar2}{Vectors of marginal GEV parameters. Required when \code{type="Qsets"} and either \code{out$method=="Bayesian"} if the marginal parameter weren't fitted or \code{"empirical"}.}
  \item{dep}{A logical value; if \code{TRUE} the estimate of the dependence is plotted when computing extreme quantile regions (\code{type="Qsets"}).}
  \item{QatCov1, QatCov2}{Matrices representing the value of the covariates at which extreme quantile regions should be computed. Required when \code{type="Qsets"}. See arguments \code{cov1} and \code{cov2} from \code{\link{fExtDep}}. }
  \item{P}{A vector indicating the probabilities associated with the quantiles to be computed. Required when \code{type="Qsets"}.}
  \item{labels}{A bivariate vector of character strings providing labels for extreme quantile regions.  Required when \code{type="Qsets"}.}
  \item{CEX}{Label and axis sizes.}
  \item{xlim, ylim}{Limits of the x and y axis when computing extreme quantile regions. Required when \code{type="Qsets"}.}
  \item{col.data, col.Qfull, col.Qfade}{Colors for data, estimate of extreme quantile regions and its credible interval (when applicable). Required when \code{type="Qsets"}.}
  \item{data}{A 2-column matrix providing the original data to be plotted when \code{type="Qsets"}.}
  \item{...}{Additional graphical parameters}
}

\details{
	If \code{type="returns"}, a (contour) plot of the probabilities of exceedances for some threshold is returned. This corresponds to the output of the \code{returns} function. \cr
	If \code{type="A"}, a plot of the estimated Pickands dependence function is drawn. If \code{A_true} is specified the plot includes the true Pickands dependence function and a functional boxplot for the estimated function.
	If \code{type="h"}, a plot of the estimated angular density function is drawn. If \code{h_true} is specified the plot includes the true angular density and a functional boxplot for the estimated function.
	If \code{type="pm"}, a plot of the prior against the posterior for the mass at \eqn{\{0\}}{\{0\}} is drawn.
	If \code{type="k"}, a plot of the prior against the posterior for the polynomial degree \eqn{k}{k} is drawn.
	If \code{type="summary"}, when the estimation was performed in a Bayesian framework then a 2 by 2 plot with types \code{"A"}, \code{"h"}, \code{"pm"} and \code{"k"} is returned. Otherwise a 1 by 2 plot with types \code{"A"} and \code{"h"} is returned.
	If \code{type="Qsets"}, extreme quantile regions are computed according to the methodology developped in Beranger et al. (2021).
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

  Marcon, G., Padoan, S.A., Naveau, P., Muliere, P., Segers, J. (2017)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials.
  \emph{Journal of Statistical Planning and Inference}, \bold{183}, 1-17.
}

\seealso{\code{\link{fExtDep.np}}.}
\examples{

###########################################################
### Example 1 - Wind Speed and Differential of pressure ###
###########################################################

data(WindSpeedGust)

years <- format(ParcayMeslay$time, format="\%Y")
attach(ParcayMeslay[which(years \%in\% c(2004:2013)),])

# Marginal quantiles
WS_th <- quantile(WS,.9)
DP_th <- quantile(DP,.9)

# Standardisation to unit Frechet (requires evd package)
pars.WS <- evd::fpot(WS, WS_th, model="pp")$estimate
pars.DP <- evd::fpot(DP, DP_th, model="pp")$estimate
  
# transform the marginal distribution to common unit Frechet:
data_uf <- trans2UFrechet(cbind(WS,DP), type="Empirical")

# compute exceedances 
rdata <- rowSums(data_uf)
r0 <- quantile(rdata, probs=.90)
extdata_WSDP <- data_uf[rdata>=r0,]

# Fit
SP_mle <- fExtDep.np(method="Frequentist", data=extdata_WSDP, k0=10, type="maxima")

# Plot
plot_ExtDep.np(out=SP_mle, type="summary")

####################################################
### Example 2 - Pollution levels in Milan, Italy ###
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
pl1 <- plot_ExtDep.np(out=pollut1, type="Qsets", summary.mcmc=pollut1_sum, 
                      mar1=gev.pars1, mar2=gev.pars2, P = 1/c(600, 1200, 2400),
                      dep=TRUE, data=data, xlim=c(0,400), ylim=c(0,400))

pl1b <- plot_ExtDep.np(out=pollut1, type="Qsets", summary.mcmc=pollut1_sum, est.out=pl1$est.out, 
                       mar1=gev.pars1, mar2=gev.pars2, P = 1/c(1200),
                       dep=FALSE, data=data, xlim=c(0,400), ylim=c(0,400))

### Frequentist estimation using Bernstein polynomials

pollut2 <- fExtDep.np(method="Frequentist", data=sdata, mar.fit=FALSE, type="rawdata", k0=8)
plot_ExtDep.np(out=pollut2, type = c("summary"), CEX=1.5)

pl2 <- plot_ExtDep.np(out=pollut2, type="Qsets", mar1=gev.pars1, mar2=gev.pars2, 
                      P = 1/c(600, 1200, 2400),
                      dep=TRUE, data=data, xlim=c(0,400), ylim=c(0,400),
                      labels=c(expression(NO[2]),expression(SO[2])),
                      col.Qfull = c("red", "green", "blue"))
                      
### Frequentist estimation using EKdH estimator

pollut3 <- fExtDep.np(method="Empirical", data=data)
plot_ExtDep.np(out=pollut3, type = c("summary"), CEX=1.5)

pl3 <- plot_ExtDep.np(out=pollut3, type="Qsets", mar1=gev.pars1, mar2=gev.pars2,
                      P = 1/c(600, 1200, 2400),
                      dep=TRUE, data=data, xlim=c(0,400), ylim=c(0,400),
                      labels=c(expression(NO[2]),expression(SO[2])),
                      col.Qfull = c("red", "green", "blue"))

}
}

\keyword{Nonparametric}
