\name{exponent_extr_mod}
\alias{exponent_extr_mod}

\title{Exponent function of extremal dependence models}
\description{Evaluate the bivariate or trivariate exponent function of the Husler-Reiss, Extremal-$t$ and Extremal Skew-$t$ models. }

\usage{
	exponent_extr_mod(model, z, param, dist)
}

\arguments{
	\item{model}{A string with the name of the model: \code{"hr"}, \code{"Extremalt"} or \code{"Skewt"}. }
	\item{z}{A vector of length \eqn{2}{2} or \eqn{3}{3}, containing strictly positive reals.}
	\item{param}{A vector containing the parameters of the model. See \bold{Details}. }	
	\item{dist}{Logical; if \code{TRUE} the distribution of the model is returned. }
}

\details{
	
	If \code{model="hr"} then the parameter vector is made of \code{choose(d,2)} positive parameters, \code{d=2,3}.
	If \code{model="Extremalt"} then the parameter vector is made of \code{choose(d,2)} dependence parameters and a degree of freedom, \code{d=2,3}.
	If \code{model="Skewt"} then the parameter vector is made of \code{choose(d,2)} dependence parameters, \code{d} shape (or skewness) parameters and a degree of freedom, \code{d=2,3}.
		
	
}

\value{
	Returns a single value correspond to the exponent function or value of the distribution.
}

\references{
	
	Beranger, B. and Padoan, S. A. (2015).
	Extreme dependence models, chapater of the book \emph{Extreme Value Modeling and Risk Analysis: Methods and Applications},
	\bold{Chapman Hall/CRC}.
	
	Beranger, B., Padoan, S. A. and Sisson, S. A. (2016).
	Models for extremal dependence derived from skew-symmetric families.
	\emph{Scandinavian Journal of Statistics}, \bold{To appear}.	
  
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

### Husler-Reiss

exponent_extr_mod(model="hr", z=c(2,3), param=1.2, dist=FALSE)
exponent_extr_mod(model="hr", z=c(2,3,1), param=c(1.2,1,1.4), dist=TRUE)

### Extremal-t

exponent_extr_mod(model="Extremalt", z=c(0.1,2), param=c(0.5,2), dist=FALSE)
exponent_extr_mod(model="Extremalt", z=c(0.1,2,3), param=c(0.5,0.4,0.9,2), dist=TRUE)

### Extremal Skew-t

exponent_extr_mod(model="Skewt", z=c(0.1,2), param=c(0.5,0,0,2), dist=FALSE)
exponent_extr_mod(model="Skewt", z=c(0.1,2,3), param=c(0.5,0.4,0.9,0,0,0,2), dist=TRUE)
exponent_extr_mod(model="Skewt", z=c(0.1,2,3), param=c(0.5,0.4,0.9,1,10,-5,2), dist=TRUE)

}




