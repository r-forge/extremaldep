\name{rExtDep}
\alias{rExtDep}

\title{Parametric and semi-parametric random generator of extreme events}
\description{This function generates random samples of iid observations from extremal dependence models and semi-parametric stochastic generators.}

\usage{
rExtDep(n, model, par, angular=FALSE, mar=c(1,1,1), num, threshold, 
        exceed.type)
}

\arguments{
	\item{n}{An integer indictaing the number of observations.}
	\item{model}{A character string with the name of the model. Parametric model include \code{"HR"} (Husler-Reiss), \code{"ET"} (Extremal-t), \code{"EST"} (Extremal Skew-t). Semi-parametric generators include \code{"semi.bvevd"} and \code{"semi.bvexceed"}. }
	\item{par}{A vector representing the parameters of the (parametric or non-parametric) model.}
	\item{angular}{A logical value; \code{TRUE} for angular outputs.}
	\item{mar}{A vector or matrix of marginal parameters.}
	\item{num}{An integer indicating the number of observations the componentwise maxima is computed over. Required when \code{model="HR"}, \code{"ET"} or \code{"EST"}. Set to \code{5e+5} unless specified.}
	\item{threshold}{A bivariate vector indicating the level of exceedances. Required when \code{model="semi.bvexceed"}.}
	\item{exceed.type}{A character string taking value "and" or "or" indicating the type of exceednaces. Required when \code{model="semi.bvexceed"}.}
}

\details{
  There is no limit of the dimensionality when \code{model="HR"}, \code{"ET"} or \code{"EST"} while \code{model="semi.bvevd"} and \code{"semi.bvexceed"} can only generate bivariate observations.
  When \code{angular=TRUE} and \code{model="semi.bvevd"} or \code{"semi.bvexceed"} the simulation of pseudo-angles follows Algorithm 1 of Marcon et al. (2017). 
  When \code{model="semi.bvevd"} and \code{angular=FALSE}, maxima samples are generated according to Algorithm 2 of Marcon et al. (2017).
  When \code{model="semi.bvexceed"} and \code{angular=FALSE}, exceedance samples are generated above value specified by \code{threshold}, according to Algorithm 3 of Marcon et al. (2017). \code{exceed.type="and"} generates samples that exceed both thresholds while \code{exceed.type="or"} generates samples exceeding at least one threshold.
  
  When the argument \code{mar} is a vector, the marginal distrutions are identical. When a matrix is provided each row corresponds to a set of marginal parameters. No marginal transformation is applied when \code{mar=c(1,1,1)}.
}

\value{
	A matrix with \eqn{n}{n} rows and \eqn{p>=2}{p>=2} columns. \eqn{p=2}{p=2} when \code{model="semi.bvevd"} or \code{"semi.bvexceed"}.   
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\references{
	Marcon, G., Naveau, P. and Padoan, S.A. (2017)
  A semi-parametric stochastic generator for bivariate extreme events
  \emph{Stat}, \bold{6}, 184-201.
}

\seealso{\code{\link{dExtDep}}, \code{\link{pExtDep}}, \code{\link{fExtDep}}, \code{\link{fExtDep.np}} }

\examples{

# Example using the trivariate Husler-Reiss
set.seed(1)
data <- rExtDep(n=10, model="HR", par=c(2,3,3))

# Example using the semi-parammetric generator of maxima
set.seed(2)
par <- ExtremalDep:::rcoef(k=6, pm=list(p0=0.05, p1=0.1))
data <- rExtDep(n=10, model="semi.bvevd", par=par$beta, 
                mar=rbind(c(0.2, 1.5, 0.6),c(-0.5, 0.4, 0.9)))

# Example using the semi-parammetric generator of maxima
set.seed(3)
data <- rExtDep(n=10, model="semi.bvexceed", par=par$beta, 
                threshold=c(0.2, 0.4), exceed.type="and")
}

\keyword{models}
