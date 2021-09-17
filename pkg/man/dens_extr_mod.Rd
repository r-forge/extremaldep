\name{dens_extr_mod}
\alias{dens_extr_mod}

\title{Density of extremal dependence models}
\description{Evaluates the bivariate density of the Extremal-$t$ and Extremal Skew-$t$ models. }

\usage{
	dens_extr_mod(model, z, param, log=TRUE)
}

\arguments{
	\item{model}{A string with the name of the model: \code{"Extremalt"} or \code{"Skewt"}. }
	\item{z}{A vector of length \eqn{2}{2}, containing strictly positive reals.}
	\item{param}{A vector containing the parameters of the model. See \bold{Details}. }
	\item{log}{Logical; if \code{TRUE} the log-density is returned. \code{FALSE}
		is the default.}
}

\details{

	If \code{model="Extremalt"} then the parameter vector is made of one dependence parameter and a degree of freedom.
	If \code{model="Skewt"} then the parameter vector is made of one dependence parameter, two shape (or skewness) parameters and a degree of freedom.

}

\value{
	Returns a single value corresponding to the density or the log-density.
}

\references{

	Beranger, B., Padoan, S. A. and Sisson, S. A. (2017).
	Models for extremal dependence derived from skew-symmetric families.
	\emph{Scandinavian Journal of Statistics}, \bold{44}(1), 21-45.

}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://mypage.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com/};
}

\examples{

### Extremal-t

dens_extr_mod(model="Extremalt", z=c(0.1,0.3), param=c(0.6,1),log=FALSE)

### Extremal Skew-t

dens_extr_mod(model="Skewt", z=c(0.1,0.3), param=c(0.6,0,0,1),log=FALSE)
dens_extr_mod(model="Skewt", z=c(0.1,0.3), param=c(0.6,-3,5,1),log=FALSE)

}

\keyword{models}
