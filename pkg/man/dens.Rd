\name{dens}
\alias{dens}

\title{Angular density and likelihood function for some extremal dependence models}
\description{Evaluates the angular density or calculates the likelihood function
	of the Pairwise Beta, Husler-Reiss, Dirichlet, Extremal-$t$, Extremal Skew-$t$ and Asymmetric Logistic models
	at one or more locations on the unit simplex. }
\usage{
dens(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), model="Pairwise", par=c(2,2,2,4), c,
        log=FALSE, vectorial=TRUE)
}
\arguments{
\item{x}{A (\eqn{n \times d}{n x d}) matrix of angular components, where the rows represent \eqn{n}{n} independent points
    in the \eqn{d}{d}-dimensional unit simplex. See \bold{Details}.
    The default is \code{rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7))}, two points in the
		\eqn{3}{3}-dimensional simplex.}
	\item{model}{A string with the name of the parametric model to be estimated. Models are
		Pairwise Beta (\code{"Pairwise"}), Husler-Reiss	(\code{"Husler"}),
		Dirichlet (\code{"Dirichlet"}), Extremal-t (\code{"Extremalt"}), Extremal Skew-t (\code{"Skewt"}) and
		Asymmetric Logistic (\code{"Asymmetric"}).}
	\item{par}{A vector containing the parameters of the model. See \bold{Details}.}
	\item{c}{A real value in \eqn{[0,1]}{[0,1]}, providing the decision rule to allocate a data point to a subset of the simplex.
		Only required for the Extremal-t, Extremal Skew-t and Asymmetric Logistic models.}
	\item{log}{Logical; if \code{TRUE} the log-density is returned. \code{FALSE}
		is the default.}
	\item{vectorial}{Logical; if \code{TRUE} when \eqn{n>1}{n>1} the different densities
		are returned as a vector of length \eqn{n}{n}. If \code{FALSE} the likelihood
		function is returned. \code{TRUE} is the default.}
}

\details{
	The Extremal-$t$ and Asymmetric Logistic models are available up to 3 dimensions;
	mass on all the subsets of the simplex is included.

	For the Pairwise Beta model, the parameter vector is decomposed as:
	\describe{
	\item{b}{A vector of size \code{choose(d,2)}. Controls the dependence between
		pairs. The default is \code{b=c(2,2,2)}.}
	\item{alpha}{A positive real that controls the general dependence between
		all the variables. The default is \eqn{4}{4}.}
	}

	For the Husler-Reiss model, the parameter vector is of size \code{choose(d,2)}.

	For the Dirichlet model, the parameter vector is decomposed a vector of size \eqn{d}{d}
	which controls the dependence between pairs.

	For the Extremal-$t$ model, the parameter vector is decomposed as:
	\describe{
	\item{rho}{A vector of size \code{choose(d,2)} representing the corrleation parameters.}
	\item{mu}{A positive integer, \eqn{\mu \geq 1}, representing the degree of freedom.}
	}

	For the Extremal Skew-$t$ model, the parameter vector is decomposed as:
	\describe{
	\item{rho}{A vector of size \code{choose(d,2)} representing the corrleation parameters.}
	\item{alpha}{A vector of size \code{d} representing the shape parameters.}
	\item{mu}{A positive integer, \eqn{\mu \geq 1}, representing the degree of freedom.}
	}

	For the Asymmetric Logistic model,  the parameter vector is decomposed as:
	\describe{
	\item{alpha}{A vector of size \eqn{1}{1} or \eqn{4}{4} depending on whether \eqn{d=2}{d=2}
		or \eqn{3}{3}.}
	\item{beta}{A vector of size \eqn{2}{2} or \eqn{9}{9} depending on whether \eqn{d=2}{d=2}
		or \eqn{3}{3}.}
	}

	If \code{log=TRUE} and \code{vectorial=FALSE} then the log-likelihood function is
	calculated.
}

\value{
	Returns a \eqn{n}{n}-dimensional vector if \code{vectorial=TRUE} or a single value
	if \code{vectorial=FALSE}.
}

\references{
	Cooley, D.,Davis, R. A., and Naveau, P. (2010).
	The pairwise beta distribution: a flexible parametric multivariate model for extremes.
	\emph{Journal of Multivariate Analysis}, \bold{101}, 2103--2117.

	Husler, J. and Reiss, R.-D. (1989),
	Maxima of normal random vectors: between independence and complete dependence,
	\emph{Statistics and Probability Letters}, \bold{7}, 283--286.

	Engelke, S., Malinowski, A., Kabluchko, Z., and Schlather, M. (2015),
	Estimation of Husler-Reiss distributions and Brown-Resnick processes,
	\emph{Journal of the Royal Statistical Society, Series B (Methodological)}, \bold{77}, 239--265.

	Coles, S. G., and Tawn, J. A. (1991),
	Modelling Extreme Multivariate Events,
	\emph{Journal of the Royal Statistical Society, Series B (Methodological)}, \bold{53}, 377--392.

	Nikoloulopoulos, A. K., Joe, H., and Li, H. (2009)
	Extreme value properties of t copulas.
	\emph{Extremes}, \bold{12}, 129--148.

	Opitz, T. (2013)
	Extremal t processes: Elliptical domain of attraction and a spectral representation.
	\emph{Jounal of Multivariate Analysis}, \bold{122}, 409--413.

	Beranger, B. and Padoan, S. A. (2015).
	Extreme dependence models, chapater of the book \emph{Extreme Value Modeling and Risk Analysis: Methods and Applications},
	\bold{Chapman Hall/CRC}.

	Beranger, B., Padoan, S. A. and Sisson, S. A. (2017).
	Models for extremal dependence derived from skew-symmetric families.
	\emph{Scandinavian Journal of Statistics}, \bold{44}(1), 21-45.

	Tawn, J. A. (1990),
	Modelling Multivariate Extreme Value Distributions,
	\emph{Biometrika}, \bold{77}, 245--253.

}


\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://mypage.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com/};
}

\examples{
if (interactive()){
### Pairwise Beta :


# Examples on the 3-dimensional simplex
# Returns the bivariate angular density at two locations

dens(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), model="Pairwise", par=c(2,2,2,4),
log=FALSE, vectorial=TRUE)

# returns the likelihood function at two locations

dens(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), model="Pairwise", par=c(2,2,2,4),
log=FALSE, vectorial=FALSE)

# returns the log-likelihood function

dens(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), model="Pairwise", par=c(2,2,2,4),
log=TRUE, vectorial=FALSE)

# Examples on the 4-dimensional simplex
# returns the bivariate angular density at two locations

dens(x=rbind(c(0.1,0.3,0.3,0.3),c(0.1,0.2,0.3,0.4)), model="Pairwise", par=c(2,2,2,1,0.5,3,4),
log=FALSE, vectorial=TRUE)

# returns the likelihood function at two locations
dens(x=rbind(c(0.1,0.3,0.3,0.3),c(0.1,0.2,0.3,0.4)), model="Pairwise", par=c(2,2,2,1,0.5,3,4),
log=FALSE, vectorial=FALSE)

# returns the log-likelihood function

dens(x=rbind(c(0.1,0.3,0.3,0.3),c(0.1,0.2,0.3,0.4)), model="Pairwise", par=c(2,2,2,1,0.5,3,4),
log=TRUE, vectorial=FALSE)


### Husler-Reiss


# Example on the 2-dimensional simplex
# returns the log-likelihood at two locations

dens(x=rbind(c(0.1,0.9),c(0.3,0.7)), model="Husler", par=1.7,
log=TRUE, vectorial=FALSE)

# Example on the 3-dimensional simplex
# returns the likelihood function at two locations

dens(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), model="Husler", par=c(1.7,0.7,1.1),
log=FALSE, vectorial=FALSE)

# Example on the 4-dimensional simplex
# returns the bivariate angular density at two locations

dens(x=rbind(c(0.1,0.1,0.4,0.4),c(0.1,0.2,0.3,0.4)), model="Husler", par=rep(1,6),
log=FALSE, vectorial=TRUE)


### Dirichlet


# Example on the 2-dimensional simplex
# returns the log-likelihood at two points

dens(x=rbind(c(0.1,0.9),c(0.3,0.7)), model="Dirichlet", par=c(1.7,0.7),
log=TRUE, vectorial=FALSE)


# Example on the 3-dimensional simplex
# returns the likelihood function at three locations

dens(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), model="Dirichlet", par=c(1.7,0.7,1.1),
log=FALSE, vectorial=FALSE)

# Example on the 4-dimensional simplex
# returns the bivariate angular density at two locations

dens(x=rbind(c(0.1,0.1,0.4,0.4),c(0.1,0.2,0.3,0.4)), model="Dirichlet", par=c(1.7,0.7,1.1,0.1),
log=FALSE, vectorial=TRUE)


### Extremal-t


# Example on the 2-dimensional simplex
# Returns the log-likelihood

dens(x=rbind(c(0.4,0.6),c(0.3,0.7)), model="Extremalt", par=c(0.7,2), c=0.1,
log=TRUE, vectorial=FALSE)

# Density in the corner

dens(x=c(0.08,0.92), model="Extremalt", par=c(0.7,2), c=0.1,
log=FALSE, vectorial=FALSE)


# Example on the 3-dimensional simplex
# Returns the log-likelihood

dens(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), model="Extremalt", par=c(rep(0.1,3),2), c=0.03,
log=FALSE, vectorial=FALSE)

# Returns the evalutaion of the angular density at three locations:
# The first one is set to be on the edge linking the second and third components
# The second one is set to be on the interior of the simplex
# The third one is set to be on the corner near the third component

if (interactive()){
dens(x=rbind(c(0.001,0.3,0.699),c(0.1,0.2,0.7),c(0.001,0.001,0.998)),
model="Extremalt", par=c(rep(0.1,3),2), c=0.01, log=FALSE, vectorial=TRUE)
}

### Extremal Skew-t


# Example on the 2-dimensional simplex
# Returns the log-likelihood

dens(x=rbind(c(0.4,0.6),c(0.3,0.7)), model="Skewt", par=c(0.7,0,0,2), c=0.1,
log=TRUE, vectorial=FALSE)

dens(x=rbind(c(0.4,0.6),c(0.3,0.7)), model="Skewt", par=c(0.7,2,-1,2), c=0.1,
log=TRUE, vectorial=FALSE)

# Density in the corner

dens(x=c(0.08,0.92), model="Skewt", par=c(0.7,0,0,2), c=0.1,
log=FALSE, vectorial=FALSE)

dens(x=c(0.08,0.92), model="Skewt", par=c(0.7,-1,2,2), c=0.1,
log=FALSE, vectorial=FALSE)

# Example on the 3-dimensional simplex
# Returns the log-likelihood

dens(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), model="Skewt", par=c(rep(0.1,3),rep(0,3),2), c=0.03,
log=FALSE, vectorial=FALSE)

# Returns the evalutaion of the angular density at three locations:
# The first one is set to be on the edge linking the second and third components
# The second one is set to be on the interior of the simplex
# The third one is set to be on the corner near the third component

if (interactive()){
dens(x=rbind(c(0.001,0.3,0.699),c(0.1,0.2,0.7),c(0.001,0.001,0.998)),
model="Skewt", par=c(rep(0.1,3),rep(0,3),2), c=0.01, log=FALSE, vectorial=TRUE)
}


### Asymmetric Logistic

# Example on the 3-dimensional simplex
# Returns the angular density at three points:
# The first one is set to be on the edge linking the second and third components
# The second one is set to be on the interior of the simplex
# The third one is set to be on the corner near the third component

dens(x=rbind(c(0.001,0.3,0.699),c(0.1,0.2,0.7),c(0.001,0.001,0.998)), c=0.05,
model="Asymmetric", par=c(1.2,1.8,4,2,rep(0.3,9)), log=FALSE, vectorial=TRUE)


}
}

\keyword{models}
