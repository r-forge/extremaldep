\name{convexity}
\alias{convexity}
\title{Convexity Constraint}

\description{
  Generate a matrix with constraints of the Bernstein polynomial
  coefficients. This allows one to obtain a polynomial representation of convex functions.}

\usage{
  convexity(v,d)}
  
\arguments{
  \item{v}{A matrix with the index set, see \bold{Details}.}
  \item{d}{The number of variables. A positive integer greater than or equal to two.}
}

\details{
  When applying the convexity constraint to the Bernstein polynomials,
  for any positive integers \code{k} and \code{d}, each row of the
  matrix (with the index set) is an ordered selection with repetition of elements in
  \eqn{0,\ldots, k}{0,...,k}. The sum over rows must be equal to
  \code{k}, that is the order of the polynomial.
  The number of rows of the index set matrix, indeed, corresponds to the
  number of coefficients \eqn{p(k,d)} involved in the Bernstein
  representation (see equation (8) in Marcon et al., 2014 for details).}

\value{
  A matrix.
}

\references{
  Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. (2014)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials. \emph{Submitted}.
}

\seealso{\code{\link{beed}}, \code{\link{bp}}, \link{beed.confband}.}

\examples{

v <- as.matrix(0:8, ncol = 1)
convexity(v, d = 2)

}


\keyword{Convexity, Bernstein polynomials}