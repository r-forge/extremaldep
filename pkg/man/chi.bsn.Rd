\name{chi.bsn}
\alias{chi.bsn}

\title{Tail dependence coefficient for the skew-normal distirbution}
\description{Evaluate the upper and lower tail dependence coefficients for the bivariate skew-normal. }

\usage{
	chi.bsn(u, corr=0, shape=rep(0,2), tail="upper")
}

\arguments{
	\item{u}{a real stricly between \eqn{0}{0} and \eqn{1}{1}.}
	\item{corr}{the correlation parameter, between \eqn{-1}{-1} and \eqn{1}{1}. }
	\item{shape}{a numeric skewness vector of length \eqn{2}{2}. }
	\item{tail}{the string \code{"upper"} or \code{"lower"}. }
}

\details{
	Approximation, the tail dependence is obtained in the limiting case where \code{u}{u} goes to eqn{1}{1}.
}

\value{
	Returns a value that is trictly greater than \eqn{0}{0} and less than \eqn{1}{1} 
	for the upper coefficient, and between \eqn{-1}{-1} and \eqn{1}{1} for the lower coefficient.
}

\references{
		
	Bortot, P.
	Tail dependence in bivariate skew-normal and skew-t distributions.
	\emph{Unpublished}.	
	
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.bocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
}

\examples{

### Upper tail dependence

chi.bsn(u=0.9,corr=0.5, shape=c(1,-2), tail="upper")

### Lower tail dependence

chi.bsn(u=0.9, corr=0.5, shape=c(1,-2), tail="lower")

}




