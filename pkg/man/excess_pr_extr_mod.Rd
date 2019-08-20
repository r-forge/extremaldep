\name{excess_pr_extr_mod}
\alias{excess_pr_extr_mod}

\title{Exceedance Probability for extremal dependence models}
\description{Exceedance Probability for bivariate or trivariate Husler-Reiss, Extremal-$t$ and Extremal Skew-$t$ models. }

\usage{
	excess_pr_extr_mod(model, z, param)
}

\arguments{
	\item{model}{A string with the name of the model: \code{"hr"}, \code{"Extremalt"} or \code{"Skewt"}. }
	\item{z}{A vector of length \eqn{2}{2} or \eqn{3}{3}, containing strictly positive reals.}
	\item{param}{A vector containing the parameters of the model. See \bold{Details}. }
}

\details{

	If \code{model="hr"} then the parameter vector is made of \code{choose(d,2)} positive parameters, \code{d=2,3}.
	If \code{model="Extremalt"} then the parameter vector is made of \code{choose(d,2)} dependence parameters and a degree of freedom, \code{d=2,3}.
	If \code{model="Skewt"} then the parameter vector is made of \code{choose(d,2)} dependence parameters, \code{d} shape (or skewness) parameters and a degree of freedom, \code{d=2,3}.


}

\value{
	Returns a probability.
}

\references{

	Beranger, B. and Padoan, S. A. (2015).
	Extreme dependence models, chapater of the book \emph{Extreme Value Modeling and Risk Analysis: Methods and Applications},
	\bold{Chapman Hall/CRC}.

	Beranger, B., Padoan, S. A. and Sisson, S. A. (2017).
	Models for extremal dependence derived from skew-symmetric families.
	\emph{Scandinavian Journal of Statistics}, \bold{44}(1), 21-45.

}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.unibocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{http://www.borisberanger.com};
}

\examples{

### Husler-Reiss

\dontrun{
excess_pr_extr_mod(model="hr", z=c(1,3), param=0.5)
excess_pr_extr_mod(model="hr", z=c(1,3,5), param=c(5,4,2))
excess_pr_extr_mod(model="hr", z=c(0.001,0.001,0.001), param=c(5,4,2))
}

### Extremal-t

\dontrun{
excess_pr_extr_mod(model="Extremalt", z=c(0.1,0.3), param=c(0.5,2))
excess_pr_extr_mod(model="Extremalt", z=c(1,3,5), param=c(0.5,0.4,0.8,2))
excess_pr_extr_mod(model="Extremalt", z=c(0.001,0.001,0.001), param=c(0.5,0.4,0.8,2))
}

### Extremal Skew-t

\dontrun{
excess_pr_extr_mod(model="Skewt", z=c(0.1,0.3), param=c(0.5,0,0,2))
excess_pr_extr_mod(model="Skewt", z=c(0.1,0.3), param=c(0.5,-10,-4,2))
excess_pr_extr_mod(model="Skewt", z=c(1,3,5), param=c(0.5,0.4,0.8,0,0,0,2))
excess_pr_extr_mod(model="Skewt", z=c(1,3,5), param=c(0.5,0.4,0.8,1,5,10,2))
excess_pr_extr_mod(model="Skewt", z=c(0.001,0.001,0.001), param=c(0.5,0.4,0.8,1,5,10,2))
}

}

\keyword{models}
