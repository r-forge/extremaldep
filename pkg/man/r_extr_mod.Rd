\name{r_extr_mod}
\alias{r_extr_mod}

\title{Random sample generation from extremal dependence models}
\description{Generates random samples of iid observations from the Extremal-$t$ and Extremal Skew-$t$ models. }

\usage{
	r_extr_mod(model, n, param)
}

\arguments{
	\item{model}{A string with the name of the model: \code{"Extremalt"} or \code{"Skewt"}. }
	\item{n}{An integer indicating the number of numbers to be generated.}
	\item{param}{A vector containing the parameters of the model. See \bold{Details}. }
}

\details{

	If \code{model="Extremalt"} then the parameter vector is made of a dependence parameter vector of size \eqn{choose(dim,2)}{choose(dim,2)} and a degree of freedom.
	If \code{model="Skewt"} then the parameter vector is made of a dependence parameter vector of size \eqn{choose(dim,2)}{choose(dim,2)}, a vector of shape (or skewness) parameters of size \eqn{dim}{dim} and a degree of freedom.

}

\value{
	Returns a matrix with \eqn{dim}{dim} columns and \eqn{n}{n} columns.
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

if (interactive()){
## Bivariate
r_extr_mod("Extremalt", n=5, par=c(0.5,2))
}

if (interactive()){
## Trivariate
r_extr_mod("Extremalt", n=5, par=c(0.5,0.6,0.4,2))
}

### Extremal Skew-t

if (interactive()){
## Bivariate
r_extr_mod("Skewt", n=5, par=c(0.5,-1,1,2))
}

if (interactive()){
## Trivariate
r_extr_mod("Skewt", n=5, par=c(0.5,0.6,0.4,-2,-2,5,2))
}

}

\keyword{models}
