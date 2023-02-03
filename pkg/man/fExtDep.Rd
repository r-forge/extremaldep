\name{fExtDep}
\alias{fExtDep}

\title{Extremal dependence estimation}
\description{This function estimates the parameters of extremal dependence models.}

\usage{
fExtDep(method="PPP", data, model, par.start = NULL, 
        c = 0, optim.method = "BFGS", trace = 0, sig = 3,
        Nsim, Nbin = 0, Hpar, MCpar, seed = NULL)
}

\arguments{
	\item{method}{A character string indicating the estimation method inlcuding \code{"PPP"}, \code{"BayesianPPP"} and \code{"Composite"}.}
	\item{data}{A matrix containing the data.}
	\item{model}{A character string with the name of the model. When \code{method="PPP"} or \code{"BayesianPPP"}, this includes \code{"PB"}, \code{"HR"}, \code{"ET"}, \code{"EST"}, \code{TD} and \code{AL} whereas when \code{method="composite"} it is restricted to \code{"HR"}, \code{"ET"} and \code{"EST"}.}
	\item{par.start}{A vector representing the initial parameters values for the optimization algorithm.}
	\item{c}{A real value in \eqn{[0,1]}{[0,1]} required when \code{method="PPP"} or \code{"BayesianPPP"} and \code{model="ET"}, \code{"EST"} and \code{"AL"}. See \code{\link{dExtDep}} for more details.}
	\item{optim.method}{A character string indicating the optimization algorithm. Required when \code{method="PPP"} or \code{"Composite"}. See \code{\link{optim}} for more details.}
	\item{trace}{A non-negative integer, tracing the progress of the optimization. Required when \code{method="PPP"} or \code{"Composite"}. See \code{\link{optim}} for more details.}
	\item{sig}{An integer indicating the number of significant digits when reporting outputs.}
	\item{Nsim}{An integer indicating the number of MCMC simulations. Required when \code{method="BayesianPPP"}.}
	\item{Nbin}{An integer indicating the length of the burn-in period. Required when \code{method="BayesianPPP"}.}
	\item{Hpar}{A list of hyper-parameters. See 'details'. Required when \code{method="BayesianPPP"}.}
	\item{MCpar}{A positive real representing the variance of the proposal distirbution. See 'details'. Required when \code{method="BayesianPPP"}.}
	\item{seed}{An integer indicating the seed to be set for reproducibility, via the routine \code{set.seed}.}
}

\details{
  When \code{method="PPP"} the approximate likelihood is used to estimate the model parameters. It relies on the \code{dExtDep} function with argument \code{method="Parametric"} and \code{angular=TRUE}.
  
  When \code{method="BayesianPPP"} a Bayesian estimation procedure of the spatral measure is considered, following Sabourin et al. (2013) and Sabourin & Naveau (2014). The argument \code{Hpar} is required to specify the hyper-parameters of the prior distributions, taking the following into consideration: 
  
  \itemize{
		\item For the \bold{Pairwise Beta model}, the parameters components are independent, log-normal.
	The vector of parameters is of size \code{choose(dim,2)+1} with positive components. The first elements are the
	pairiwse dependence parameters \code{b} and the last one is the global dependence parameter \code{alpha}.
	The list of hyper-parameters should be of the form
	\code{mean.alpha=, mean.beta=, sd.alpha=, sd.beta=};
		\item For the \bold{Husler-Reiss model}, the parameters are independent, log-normally distributed.
	The elements correspond to the \code{lambda} parameter. The list of hyper-parameters should be of the form \code{mean.lambda=, sd.lambda=};
		\item For the \bold{Dirichlet model}, the parameters are independent, log-normally distributed.
	The elements correspond to the \code{alpha} parameter. The list of hyper-parameters should be of the form \code{mean.alpha=, sd.alpha=};
		\item For the \bold{Extremal-t model}, the parameters are independent, logit-squared for \code{rho} and log-normal for \code{mu}. The first elements correspond	to the correlation parameters \code{rho} and the last parameter is the global dependence parameter \code{mu}. The list of hyper-parameters should be of the form \code{mean.rho=, mean.mu=, sd.rho=, sd.mu=};
		\item For the \bold{Extremal skewt-t model}, the parameters are independent, logit-squared for \code{rho}, normal for \code{alpha} and log-normal for \code{mu}. The first elements correspond	to the correlation parameters \code{rho}, then the skewness parameters \code{alpha} and the last parameter is the global dependence parameter \code{mu}. The list of hyper-parameters should be of the form \code{mean.rho=, mean.alpha=, mean.mu=, sd.rho=, sd.alpha=, sd.mu=};
		\item For the \bold{Asymmetric Logistic model}, the parameters' components are independent, log-normal for \code{alpha} and logit for \code{beta}. The list of hyper-parameters should be of the form \code{mean.alpha=, mean.beta=, sd.alpha=, sd.beta=}.
	}
	
	The proposal distribution for each (transformed) parameter is a normal distribution centred on the (transformed) current parameter value, with variance \code{MCpar}. 
	
	When \code{method="Composite"}, the pairwise composite likelihood is applied, based on the \code{dExtDep} function with argument \code{method="Parametric"} and \code{angular=FALSE}.
  
}

\value{

  When \code{method == "PPP"} or \code{"Composite"}, a list is returned including
  \itemize{
    \item{\code{par}: }{The estimated parameters.}
    \item{\code{LL}: }{The maximised log-likelihood.}
    \item{\code{SE}: }{The standard errors.}
    \item{\code{TIC}: }{The Takeuchi Information Criterion.}
  }
  
  When \code{method == "BayesianPPP"}, a list is returned including 
  
  \itemize{
  	\item{\code{stored.vales}: }{A \eqn{(Nsim-Nbin)*d}{(Nsim-Nbin)*d} matrix, where \eqn{d}{d} is the dimension of the parameter space}
  	\item{\code{llh}: }{A vector of size \eqn{(Nsim-Nbin)}{(Nsim-Nbin)} containing the log-likelihoods evaluadted at each parameter
		of the posterior sample.}
	  \item{\code{lprior}: }{A vector of size \eqn{(Nsim-Nbin)}{(Nsim-Nbin)} containing the logarithm of the prior densities evaluated
		at each parameter of the posterior sample.}
		\item{\code{arguments}: }{The specifics of the algorithm.}
	  \item{\code{elapsed}: }{The time elapsed, as given by \code{proc.time} between the start and end of the run.}
	  \item{\code{Nsim}: }{The same as the passed argument.}
	  \item{\code{Nbin}: }{Idem.}
	  \item{\code{n.accept}: }{The total number of accepted proposals.}
	  \item{\code{n.accept.kept}: }{The number of accepted proposals after the burn-in period. }
	  \item{\code{emp.mean}: }{The estimated posterior parameters mean.}
	  \item{\code{emp.sd}: }{The empirical posterior sample standard deviation.}
	  \item{\code{BIC}: }{The Bayesian Information Criteria.}
  }
              
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\references{
	Beranger, B. and Padoan, S. A. (2015).
	Extreme dependence models, chapater of the book \emph{Extreme Value Modeling and Risk Analysis: Methods and Applications},
  \bold{Chapman Hall/CRC}.
  
  Sabourin, A., Naveau, P. and Fougeres, A-L (2013)
  Bayesian model averaging for multivariate extremes
  \emph{Extremes}, \bold{16}, 325-350.
  
  Sabourin, A. and Naveau, P. (2014)
  Bayesian Dirichlet mixture model for multivariate extremes: A re-parametrization
  \emph{Computational Statistics & Data Analysis}, \bold{71}, 542-567.
}

\seealso{\code{\link{dExtDep}}, \code{\link{pExtDep}}, \code{\link{rExtDep}}, \code{\link{fExtDep.np}} }

\examples{

# Example using the Poisson Point Proce Process appraoch
data(pollution)
if (interactive()){
  f.hr <- fExtDep(method="PPP", data=PNS, model="HR", 
                  par.start = rep(0.5, 3), trace=2)
}

# Example using the pairwise composite (full) likelihood
if (interactive()){
  set.seed(1)
  data <- rExtDep(n=300, model="ET", par=c(0.6,3))
  f.et <- fExtDep(method="Composite", data=data, model="ET", 
                  par.start = c(0.5, 1), trace=2)
}
}

\keyword{models}
