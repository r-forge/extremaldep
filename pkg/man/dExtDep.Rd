\name{dExtDep}
\alias{dExtDep}

\title{Parametric and non-parametric density of Extremal Dependence}
\description{This function calculates the density of parametric multivariate extreme distributions and corresponding angular density, or the non-parametric angular density represented through Bernstein polynomials.}

\usage{
dExtDep(x, method="Parametric", model, par, angular=TRUE, log=FALSE, 
        c=NULL, vectorial=TRUE, mixture=FALSE)
}

\arguments{
	\item{x}{A vector or a matrix. The value at which the density is evaluated.}
	\item{method}{A character string taking value \code{"Parametric"} or \code{"NonParametric"}}
	\item{model}{A string with the name of the model: \code{"PB"} (Pairwise Beta), \code{"HR"} (Husler-Reiss), \code{"ET"} (Extremal-t), \code{"EST"} (Extremal Skew-t), \code{TD} (Tilted Dirichlet) or \code{AL} (Asymmetric Logistic) when evaluating the angular density. Restricted to \code{"HR"}, \code{"ET"} and \code{"EST"} for multivariate extreme value densities. Required when \code{method="Parametric"}.}
	\item{par}{A vector representing the parameters of the (parametric or non-parametric) model.}
	\item{angular}{A logical value specifying if the angular density is computed.}
	\item{log}{A logical value specifying if the log density is computed.}
	\item{c}{A real value in \eqn{[0,1]}{[0,1]}, providing the decision rule to allocate a data point to a subset of the simplex.
		Only required for the parametric angular density of the Extremal-t, Extremal Skew-t and Asymmetric Logistic models.}
	\item{vectorial}{A logical value; if \code{TRUE} a vector of \code{nrow(x)} densities is returned. If \code{FALSE} the likelihood
		function is returned (product of densities or sum if \code{log=TRUE}. \code{TRUE} is the default.}
	\item{mixture}{A logical value specifying if the Bernstein polynomial representation of distribution should be expressed as a mixture. If \code{mixture=TRUE} the density integrates to 1. Required when \code{method="NonParametric"}. }
		
}

\details{
  Note that when \code{method="Parametric"} and \code{angular=FALSE}, the density is only available in 2 dimensions.
  When \code{method="Parametric"} and \code{angular=TRUE}, the models \code{"AL"}, \code{"ET"} and \code{"EST"} are limited to 3 dimensions. This is because of the existence of mass on the subspaces on the simplex (and therefore the need to specify \code{c}).

  For the parametric models, the appropriate length of the parameter vector can be obtained from the \code{dim_ExtDep} function and are summarized as follows.
  When \code{model="HR"}, the parameter vector is of length \code{choose(dim,2)}.
  When \code{model="PB"} or \code{model="Extremalt"}, the parameter vector is of length \code{choose(dim,2) + 1}.
  When \code{model="EST"}, the parameter vector is of length \code{choose(dim,2) + dim + 1}.
  When \code{model="TD"}, the parameter vector is of length \code{dim}.
  When \code{model="AL"}, the parameter vector is of length \code{2^(dim-1)*(dim+2) - (2*dim+1)}.
}

\value{
	If \code{x} is a matrix and \code{vectorial=TRUE}, a vector of length \code{nrow(x)}, otherwise a scalar.   
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\references{

	Beranger, B. and Padoan, S. A. (2015).
	Extreme dependence models, chapater of the book \emph{Extreme Value Modeling and Risk Analysis: Methods and Applications},
	\bold{Chapman Hall/CRC}.

	Beranger, B., Padoan, S. A. and Sisson, S. A. (2017).
	Models for extremal dependence derived from skew-symmetric families.
	\emph{Scandinavian Journal of Statistics}, \bold{44}(1), 21-45.

	Coles, S. G., and Tawn, J. A. (1991),
	Modelling Extreme Multivariate Events,
	\emph{Journal of the Royal Statistical Society, Series B (Methodological)}, \bold{53}, 377--392.

	Cooley, D.,Davis, R. A., and Naveau, P. (2010).
	The pairwise beta distribution: a flexible parametric multivariate model for extremes.
	\emph{Journal of Multivariate Analysis}, \bold{101}, 2103--2117.

	Engelke, S., Malinowski, A., Kabluchko, Z., and Schlather, M. (2015),
	Estimation of Husler-Reiss distributions and Brown-Resnick processes,
	\emph{Journal of the Royal Statistical Society, Series B (Methodological)}, \bold{77}, 239--265.

	Husler, J. and Reiss, R.-D. (1989),
	Maxima of normal random vectors: between independence and complete dependence,
	\emph{Statistics and Probability Letters}, \bold{7}, 283--286.

	Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. and Segers, J. (2017)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials.
  \emph{Journal of Statistical Planning and Inference}, \bold{183}, 1-17.

	Nikoloulopoulos, A. K., Joe, H., and Li, H. (2009)
	Extreme value properties of t copulas.
	\emph{Extremes}, \bold{12}, 129--148.

	Opitz, T. (2013)
	Extremal t processes: Elliptical domain of attraction and a spectral representation.
	\emph{Jounal of Multivariate Analysis}, \bold{122}, 409--413.

	Tawn, J. A. (1990),
	Modelling Multivariate Extreme Value Distributions,
	\emph{Biometrika}, \bold{77}, 245--253.

}

\seealso{\code{\link{pExtDep}}, \code{\link{rExtDep}}, \code{\link{fExtDep}}, \code{\link{fExtDep.np}} }

\examples{

# Example of PB on the 4-dimensional simplex
dExtDep(x=rbind(c(0.1,0.3,0.3,0.3),c(0.1,0.2,0.3,0.4)), method="Parametric", 
        model="PB", par=c(2,2,2,1,0.5,3,4), log=FALSE)

# Example of EST in 2 dimensions
dExtDep(x=c(1.2,2.3), method="Parametric", model="EST", par=c(0.6,1,2,3), angular=FALSE, log=TRUE)

# Example of non-parametric angular density
par <- ExtremalDep:::rcoef(k=6, pm=list(p0=0.05, p1=0.1))
dExtDep(x=rbind(c(0.1,0.9),c(0.2,0.8)), method="NonParametric", par=par$beta)

}

\keyword{models}
