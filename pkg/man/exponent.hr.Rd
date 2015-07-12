\name{exponent.hr}
\alias{exponent.hr}
\alias{print.exponent.hr}

\title{Exponent function for parametric models.}

\description{Exponent function of selected model.}

\usage{
	exponent.hr(X, par, dist=FALSE)
}

\arguments{
	\item{X}{The values, on Frechet scale, at which the exponent function is evaluated. 
		The length of the vector indicates the dimension (dim) considered.}
	\item{par}{A vector of size \code{choose(dim,2)} representing the parameters of the model.}
	\item{dist}{Logical; if \code{dist=TRUE} then the distribution function of the model is returned.
		Deault is \code{dist=FALSE} which return the exponent function.}
}

\details{The maximum dimension is \eqn{5}{5}. See Beranger and Padoan (2015) for the full expression.
	If X is a vector of ones then we obtain the coefficient of tail dependence.}
  
\references{  
	Beranger, B. and Padoan, S. A. (2015).
	Extreme dependence models, chapater of the book \emph{Extreme Value Modeling and Risk Analysis: Methods and Applications},
  \bold{Chapman Hall/CRC}.
}