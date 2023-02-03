\name{fGEV}
\alias{fGEV}

\title{Fitting of the Generalized Extreme Value Distribution}
\description{Maximum-likelihood and Metropolis-Hastings algorithm for the estimation of the generalized extreme value distribution.}

\usage{
fGEV(data, par.start, method="Frequentist", u, cov, 
     optim.method="BFGS", optim.trace=0, sig0, nsim)
}

\arguments{
	\item{data}{A vector representing the data, which may contain missing values.}
	\item{par.start}{A vector of length \eqn{3}{3} giving the starting values for the parameters to be estimated. If missing, moment estimates are used.}
	\item{method}{A character string indicating whether the estimation is done following a \code{"Frequentist"} or \code{"Bayesian"} approach.}
	\item{u}{A real indicating a high threshold. If supplied a threshld exceedance approach is taken and computations use the censored likelihood. If missing, a block maxima approach is taken and the regular GEV likelihood is used.}
	\item{cov}{A matrix of covariates to define a linear model for the location parameter.}
	\item{optim.method}{The optimization method to be used. Required when \code{method="Frequentist"}. See \code{\link{optim}} for more details.}
	\item{optim.trace}{A non-negative interger tracing the progress of the optimization. Required when \code{method="Frequentist"}. See \code{\link{optim}} for more details.}
	\item{sig0}{Positive reals representing the initial value for the scaling parameter of the multivariate normal proposal distribution for both margins. Required when \code{method="Bayesian"}.}
	\item{nsim}{An integer indicating the number of iterations in the Metropolis-Hastings algorithm. Required when \code{method="Bayesian"}.}
}

\details{
  When \code{cov} is a vector of ones then the location parameter \eqn{\mu}{\mu} is constant. On the contrary, when \code{cov} is provided, it represents the design matrix for the linear model on \eqn{\mu}{\mu} (the number of columns in the matrix indicates the number of linear predictors).
  When \code{u=NULL} or missing, the likelihood function is given by \deqn{\prod_{i=1}^{n}g(x_i; \mu, \sigma, \xi)}{\prod_{i=1}^{n}f(x_i; \mu, \sigma, \xi)} where \eqn{g(\cdot;\mu,\sigma,\xi)}{g(\cdot;\mu,\sigma,\xi)} represents the GEV pdf, whereas when a threshold value is set the likelihood is given by \deqn{k_n \log\left( G(u;\mu,\sigma,\xi) \right) \times \prod_{i=1}^n \frac{\partial}{\partial x}G(x;\mu,\sigma,\xi)|_{x=x_i} }{k_n \log\left( G(u;\mu,\sigma,\xi) \right) \times \prod_{i=1}^n \frac{\partial}{\partial x}G(x;\mu,\sigma,\xi)|_{x=x_i} } where \eqn{G(\cdot;\mu,\sigma,\xi)}{G(\cdot;\mu,\sigma,\xi)} is the GEV cdf and \eqn{k_n}{k_n} is the empirical estimate of the probability of being greater than the threshold \code{u}. 
  
  Note that the case \eqn{\xi<=0}{\xi<=0} is not yet considered when \code{u} is considered.
  
  The choice \code{method="Bayesian"} applies a random walk Metropolis-Hastings algorithm as described in Section 3.1 and Step 1 and 2 of Algorithm 1 from Beranger et al. (2021). The algorithm may restart for several reasons including if the proposed value of the parameters changes too much from the current value (see Garthwaite et al. (2016) for more details.)
  
  The choice \code{method="Frequentist"} uses the \code{\link{optim}} function to find the maximum likelihood estimator.
}

\value{
	When \code{method="Frequentist"} the routine returns a list including the parameter estimates (\code{est}) and associated variance-covariance matrix (\code{varcov}) and standard errors (\code{stderr}).
	When \code{method="Bayesian"} the routine returns a list including
	\itemize{
	\item{\code{param_post}}{the parameter posterior sample;}
	\item{\code{accepted}}{a binary vector indicating which proposal was accepted;} 
	\item{\code{straight.reject}}{a binary vector indicating which proposal were rejected straight away given that the proposal is a multivariate normal and there are conditions on the parameter values;}
	\item{\code{nsim}}{the number of simulations in the algorithm;}
	\item{\code{sig.vec}}{the vector of updated scaling parameter in the multivariate normal porposal distribution at each iteration;}
	\item{\code{sig.restart}}{the value of the scaling parameter in the multivariate normal porposal distribution when the algorithm needs to restart;}
	\item{\code{acc.vec}}{a vector of acceptance probabilities at each iteration.}
	}
	
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\references{
	Beranger, B., Padoan, S. A. and Sisson, S. A. (2021).
	Estimation and uncertainty quantification for extreme quantile regions.
  \emph{Extremes}, \bold{24}, 349-375.
  
  Garthwaite, P. H., Fan, Y. and Sisson S. A. (2016). 
  Adaptive optimal scaling of Metropolis-Hastings algorithms using the Robbins-Monro process. 
  \emph{Communications in Statistics - Theory and Methods}, \bold{45}(17), 5098-5111.
}  

\seealso{ \code{\link{dGEV}} }
  
\examples{

##################################################
### Example - Pollution levels in Milan, Italy ###
##################################################

data(MilanPollution)

# Frequentist estimation
fit <- fGEV(Milan.winter$PM10)
fit$est

# Bayesian estimation with high threshold
cov <- cbind(rep(1,nrow(Milan.winter)), Milan.winter$MaxTemp, 
             Milan.winter$MaxTemp^2)
u <- quantile(Milan.winter$PM10, prob=0.9, type=3, na.rm=TRUE)
if(interactive()){
  fit2 <- fGEV(data=Milan.winter$PM10, par.start=c(50,0,0,20,1), 
               method="Bayesian", u=u, cov=cov, sig0=0.1, nsim=5e+4) 
}
}


