\name{dim_ExtDep}
\alias{dim_ExtDep}

\title{Dimensions calculations for parametric extremal dependence models}
\description{This function calculates the dimensions of an extremal dependence model for a given set of parameters, the dimension of the parameter vector for a given dimension and verifies the adequacy between model dimension and length of parameter vector when both are provided.}

\usage{
dim_ExtDep(model, par=NULL, dim=NULL)
}

\arguments{
	\item{model}{A string with the name of the model: \code{"PB"} (Pairwise Beta), \code{"HR"} (Husler-Reiss), \code{"ET"} (Extremal-t), \code{"EST"} (Extremal Skew-t), \code{TD} (Tilted Dirichlet) or \code{AL} (Asymmetric Logistic). }
	\item{par}{A vector representing the parameters of the model.}
	\item{dim}{An integer representing the dimension of the model. }
}

\details{

  One of \code{par} or \code{dim} need to be provided.
  If \code{par} is provided, the dimension of the model is calculated.
  If \code{dim} is provided, the length of the parameter vector is calculated.
  If both \code{par} and \code{dim} are provided, the adequacy between the dimension of the model and the length of the parameter vector is checked.

  For \code{model="HR"}, the parameter vector is of length \code{choose(dim,2)}.
  For \code{model="PB"} or \code{model="Extremalt"}, the parameter vector is of length \code{choose(dim,2) + 1}.
  For \code{model="EST"}, the parameter vector is of length \code{choose(dim,2) + dim + 1}.
  For \code{model="TD"}, the parameter vector is of length \code{dim}.
  For \code{model="AL"}, the parameter vector is of length \code{2^(dim-1)*(dim+2) - (2*dim+1)}.
}

\value{
	If \code{par} is not provided and \code{dim} is provided: returns an integer indicating the length of the parameter vector.
	If \code{par} is provided and \code{dim} is not provided: returns an integer indicating the dimension of the model.
  If \code{par} and \code{dim} are provided: returns a \code{TRUE/FALSE}statement indicating whether the length of the parameter and the dimension match.
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\examples{
  dim_ExtDep(model="EST", dim=3)
  dim_ExtDep(model="AL", dim=3)
  
  dim_ExtDep(model="PB", par=rep(0.5,choose(4,2)+1) )
  dim_ExtDep(model="TD", par=rep(1,5) )
  
  dim_ExtDep(model="EST", dim=2, par=c(0.5,1,1,1) )
  dim_ExtDep(model="PB", dim=4, par=rep(0.5,choose(4,2)+1) )
}

\keyword{models}
