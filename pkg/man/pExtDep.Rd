\name{pExtDep}
\alias{pExtDep}

\title{Parametric and non-parametric distribution function of Extremal Dependence}
\description{This function evaluates the distribution function of parametric multivariate extreme distributions and the angular probability distribution represented through Bernstein polynomials.}

\usage{
pExtDep(q, type, method="Parametric", model, par, plot=TRUE, 
        main, xlab, cex.lab, cex.axis, lwd,...)
}

\arguments{
	\item{q}{A vector or matrix of quantiles.}
	\item{type}{A character string taking value \code{"lower"}, \code{"inv.lower"} or \code{"upper"}. Required when \code{method="Parametric"}. }
	\item{method}{A character string taking value \code{"Parametric"} or \code{"NonParametric"}}
	\item{model}{A character string with the name of the model: \code{"HR"} (Husler-Reiss), \code{"ET"} (Extremal-t) or \code{"EST"} (Extremal Skew-t). Required when \code{method="Parametric"}. }
	\item{par}{A vector or a matrix representing the parameters of the (parametric or non-parametric) model. When in matrix format, rows indicate different sets of parameter values.}
	\item{plot}{A logical value; if \code{TRUE} (default) a plot is displayed. See \bold{details}.}
	\item{main, xlab, cex.lab, cex.axis, lwd}{Arguments of the \code{hist()} function.}
	\item{...}{Additional graphical parameter when \code{plot=TRUE}.}
}

\details{
  Note that when \code{method="Parametric"}, the distribution function is only available in 2 and 3 dimensions. Refer to the \code{\link{dim_ExtDep}} function for the appropriate length of the parameter vector. \cr
  When \code{type="lower"}, the cumulative distribution function is computed, i.e., \deqn{G(\boldsymbol{x}) = P\left(\boldsymbol{X} \leq \boldsymbol{x}\right), x \in R^d, d=2,3.}{G(x) = P(X<=x), x \in R^d, d=2,3.} \cr
  When \code{type="inv.lower"}, the survival function is computed, i.e., \deqn{1-G(\boldsymbol{x}) = 1-P\left(\boldsymbol{X} \leq \boldsymbol{x}\right).}{1-G(x) = 1-P(X<=x).} This corresponds to the probability of at least one component of \eqn{X}{X} is greater than its corresponding element in \eqn{x}{x}. \cr
  When \code{type="upper"}, the joint probability of exceedance is computed, i.e., \deqn{P\left(\boldsymbol{X} \geq \boldsymbol{x}\right).}{P(X>x).} \cr
  Finally, when \code{method="NonParametric"}, the distribution function is only available in 2 dimensions.\cr
  
  The argument \code{plot} is only applicable when \code{par} is a matrix. Typically its main use should be when \code{par} corresponds to some posterior sample (e.g. from \code{fExtDep} with \code{moethod="BayesianPPP"}). A histogram of the probabilities evaluated at each set of parameters is displayed, as well as a kernel density estimator, \eqn{2.5\%, 50\%, 97.5\%}{2.5\%, 50\%, 97.5\%} quantiles (crosses) and mean (dot). The argument \code{...} is used to specify additional parameters in the \code{hist()} function.
}

\value{
	When \code{par} is a vector: if \code{q} is a matrix the function returns a vector of length \code{nrow(q)}, otherwise a scalar. \cr
	When \code{par} is a matrix: if \code{q} is a matrix the function returns a matrix with \code{nrow(par)} rows and \code{nrow(q)} columns, otherwise a vector of length \code{nrow(par)}. 
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

	Husler, J. and Reiss, R.-D. (1989),
	Maxima of normal random vectors: between independence and complete dependence,
	\emph{Statistics and Probability Letters}, \bold{7}, 283--286.

	Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. and Segers, J. (2017)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials.
  \emph{Journal of Statistical Planning and Inference}, \bold{183}, 1-17.
}

\seealso{\code{\link{dExtDep}}, \code{\link{rExtDep}}, \code{\link{fExtDep}}, \code{\link{fExtDep.np}} }

\examples{

# Example using the trivariate Extremal Skew-t
pExtDep(q=c(1,1.2, 0.6), type="lower", method="Parametric", 
        model="EST", par=c(0.2, 0.4, 0.6,2,2,2,1))

# Example using the bivariate Extremal-t
pExtDep(q=rbind(c(1.2, 0.6), c(1.1, 1.3)), type="inv.lower", 
        method="Parametric", model="ET", par=c(0.2, 1))
pExtDep(q=rbind(c(1.2, 0.6), c(1.1, 1.3)), type="inv.lower", 
        method="Parametric", model="EST", par=c(0.2, 0, 0, 1))

# Example of non-parametric angular density
par <- ExtremalDep:::rcoef(k=6, pm=list(p0=0.05, p1=0.1))
pExtDep(q=rbind(c(0.1,0.9),c(0.2,0.8)), method="NonParametric", par=par$beta)

}

\keyword{models}
