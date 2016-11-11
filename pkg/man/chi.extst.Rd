\name{chi.extst}
\alias{chi.extst}

\title{Tail dependence coefficient for the Extremal Skew-$t$ model}
\description{Evaluate the upper and lower tail dependence coefficients for the bivariate Extremal Skew-$t$ model. }

\usage{
	chi.extst(corr=0, shape=rep(0,2), df=1, tail="upper")
}

\arguments{
	\item{corr}{the correlation parameter, between \eqn{-1}{-1} and \eqn{1}{1}. }
	\item{shape}{a numeric skewness vector of length \eqn{2}{2}. }
	\item{df}{a single positive value representing the degree of freedom. }
	\item{tail}{the string \code{"upper"} or \code{"lower"}. }
}


\value{
	Returns a value that is trictly greater than \eqn{0}{0} and less than \eqn{1}{1}.
}


\author{
  Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{faculty.bocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
}

\references{
		
	Padoan, S. A. (2011).
	Multivariate extreme models based on underlying skew-t and skew-normal distributions.
	\emph{Journal of Multivariate Analysis}, \bold{102}, 977-991.	
	
}

\examples{

### Upper tail dependence

chi.extst(corr=0.5, shape=c(1,-2), df=2, tail="upper")

### Lower tail dependence

chi.extst(corr=0.5, shape=c(1,-2), df=2, tail="lower")

}




