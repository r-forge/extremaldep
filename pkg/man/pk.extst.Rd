\name{pk.extst}
\alias{pk.extst}

\title{Pickands dependence function for the Extremal Skew-$t$ model.}
\description{Evaluate the bivariate and trivariate Pickands dependence function for the bivariate skew-normal. }

\usage{
	pk.extst(x, param)
}

\arguments{
	\item{x}{a vector of length \eqn{2}{2} or \eqn{3}{3} that belongs to the corresponding simplex. }
	\item{param}{the parameter vector, containing the dependence, shape and df parameters. }
}

\details{
	In the bivariate case, there is \eqn{1}{1} dependence parameter, \eqn{2}{2} shape parameters and a degree of freedom.
	In the trivariate case, there is \eqn{3}{3} dependence parameter, \eqn{3}{3} shape parameters and a degree of freedom.
	Dependence parameters must be between \eqn{-1}{-1} and \eqn{1}{1}, the degree of freedom must be positive.
}

\references{
  Padoan, S. A. (2011).
	Multivariate extreme models based on underlying skew-t and skew-normal distributions.
	\emph{Journal of Multivariate Analysis}, \bold{102}, 977-991.	
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{faculty.bocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
}

\examples{

### Upper tail dependence

pk.extst(x=c(0.5,0.5), param=c(0.4,-2,4,3))

### Lower tail dependence

pk.extst(x=c(0.2,0.4,0.4), param=c(0.4,0.3,0.7,3,-1,0,2))

}




