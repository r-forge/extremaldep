\name{excess.dist.hr}
\alias{excess.dist.hr}
\alias{print.excess.dist.hr}

\title{Tail function for parametric models.}

\description{Tail function of selected model.}

\usage{
excess.dist.hr(X, par)
}

\arguments{
	\item{X}{The values, on Frechet scale, at which the tail function is evaluated. 
		The length of the vector indicates the dimension (dim) considered.}
	\item{par}{A vector of size \code{choose(dim,2)} representing the parameters of the model}
}

\details{The maximum dimension is \eqn{5}{5}. See Beranger and Padoan (2015) for the full expression.
	If X is a vector of ones then the upper tail dependence function is calculated.
}

\author{
  Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.bocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
}

\references{	
	Beranger, B. and Padoan, S. A. (2015).
	Extreme dependence models, chapater of the book \emph{Extreme Value Modeling and Risk Analysis: Methods and Applications},
  \bold{Chapman Hall/CRC}.
}

\examples{
## Coefficient of upper tail dependence

## Load datasets
# data(pollution) 

# est.hr.PNNS <- posteriorMCMC(Nsim=50e+4,Nbin=30e+4,dat= PNNS, model='Husler',Hpar=list(mean.lambda=0, sd.lambda=3), MCpar=0.35,seed=14342)
# excess.dist.hr(rep(1,4), est.hr.PNNS$emp.mean)
}
