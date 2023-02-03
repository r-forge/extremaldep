\name{index.ExtDep}
\alias{index.ExtDep}

\title{Index of extremal dependence}
\description{This function computes the extremal coefficient, Pickands dependence function and the coefficients of upper and lower tail dependence for several parametric models and the lower tail dependence function for the bivairate skew-normal distribution.}

\usage{
index.ExtDep(object, model, par, x, u)
}

\arguments{
  \item{object}{A character string indicating the index of extremal dependence to compute, including the extremal coefficient \code{"extremal"}, the Pickands dependence function \code{"pickands"}, the coefficient of upper tail dependence \code{"upper.tail"} and the coefficient of lower tail dependence \code{"lower.tail"}.}
	\item{model}{A character string indicating the model/distribution. When \code{object="extremal"}, \code{"pickands"} or \code{"upper.tail"} corresponding quantities can be calculated for the Husler-Reiss (\code{"HR"}), extremal-t (\code{"ET"}) and extremal skew-t (\code{"EST"}) are available. When \code{object="lower.tail"} then the extremal-t (\code{"ET"}) and extremal skew-t (\code{"EST"}) models are available as well as the skew-normal distribution (\code{"SN"}).}
  \item{par}{A vector indicating the parameter values of the corresponding model/distribution.}
  \item{x}{A vector on the bivariate or trivariate unit simplex indicating where to evaluate the Pickands dependence function.}
  \item{u}{A real in \eqn{[0,1]}{[0,1]} indicating the value at which to evaluate the lower tail dependence function of the bivariate skew-normal distribution.}
}

\details{
  The extremal coefficient is defined as \deqn{\theta = V(1,\ldots,1) = d \int_{W} \max_{j \in \left\{1, ..., d\right\} }(w_j) dH(w) = - \log G(1,\ldots,1),}{\theta = V(1,...,1) = d \int_{W} \max_{j \in {1, ..., d} }(w_j) dH(w) = - \log G(1,...,1),} where \eqn{W}{W} represents the unit simplex, \eqn{V}{V} is the exponent function and \eqn{G(\cdot)}{G(.)} the distribution function of a multivariate extreme value model. Bivariate and trivariate versions are available.
  
  The Pickands dependence function is defined as \eqn{A(x) = - \log G(1/x)}{A(x) = - \log G(1/x)} for \eqn{x}{x} in the bivariate/trivariate simplex (\eqn{W}{W}). 
  
  The coefficient of upper tail dependence is defined as \deqn{\vartheta = R(1,\ldots,1) = d \int_{W} \min_{j \in \left\{1, ..., d\right\} }(w_j) dH(w).}{\vartheta = R(1,...,1) = d \int_{W} \min_{j \in {1, ..., d} }(w_j) dH(w).} In the bivariate case, using the inclusion-exclusion principle this reduces to \eqn{\vartheta =  2 + \log G(1,1) = 2 - V(1,1)}{\vartheta = 2 + \log G(1,1) = 2 - V(1,1)}.
  
  For the skew-normal distribution, the lower tail dependence function is defined as in Bortot (2010). This is an approximation where the tail dependence is obtained in the limiting case where \code{u}{u} goes to eqn{1}{1}. The \code{par} argument should be a vector of length \eqn{3}{3} comprising of the correlation parameter, between \eqn{-1}{-1} and \eqn{1}{1} and two real-valued skewness parameters.
  
}

\value{
  When \code{object="extremal"}, returns a value between \eqn{1}{1} and \eqn{d}{d} (\eqn{d=2,3}{d=2,3}). \cr
  When \code{object="pickands"}, returns a value between \eqn{\max(x)}{\max(x)} and \eqn{1}{1}. \cr
  When \code{object="upper.tail"}, returns a value between \eqn{0}{0} and \eqn{1}{1}. \cr
  When \code{object="lower.tail"}, returns a value between \eqn{-1}{-1} and \eqn{1}{1}.
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\references{
	Bortot, P. (2010)
	Tail dependence in bivariate skew-normal and skew-t distributions.
	\emph{Unpublished}.
}

\examples{

#############################
### Extremal skew-t model ###
#############################

### Extremal coefficient
index.ExtDep(object="extremal", model="EST", par=c(0.5,1,-2,2))

### Pickands dependence function
w <- seq(0.00001, .99999, length=100)
pick <- vector(length=100)
for(i in 1:100){
  pick[i] <- index.ExtDep(object="pickands", model="EST", par=c(0.5,1,-2,2), 
                          x=c(w[i],1-w[i]))
}

plot(w, pick, type="l", ylim=c(0.5, 1), ylab="A(t)", xlab="t")
polygon(c(0, 0.5, 1), c(1, 0.5, 1), lwd=2, border = 'grey')

### Upper tail dependence coefficient
index.ExtDep(object="upper.tail", model="EST", par=c(0.5,1,-2,2))

### Lower tail dependence coefficient
index.ExtDep(object="lower.tail", model="EST", par=c(0.5,1,-2,2))

################################
### Skew-normal distribution ###
################################

### Lower tail dependence function
index.ExtDep(object="lower.tail", model="SN", par=c(0.5,1,-2), u=0.5)

}

\keyword{models}
