\name{Convexity}
\alias{Convexity}
\title{Convexity Constraint}
\description{
  Generate the matrix of constraints for the \code{beta} coefficients which guarantees that the Bernstein representation of a function is convex.
}
\usage{
  Convexity(v,d)
}
\arguments{
  \item{v}{Matrix of index set (see \bold{Details}).}
  \item{d}{The dimension. An integer greater than or equal to two.}
}

\details{
  When applying the convexity constraint to the Bernstein polynomials, for any \code{k} and {d}, each row of the index set is an ordered selection with repetition of elements in \eqn{(0,\ldots, k)}{(0,...,k)}. The sum over rows must be equal to \code{k}.
  The number of rows of the index set matrix, indeed, corresponds to the number of coefficients \eqn{p(k,d)} involved in the Bernstein representation (the relation is explicated in equation (8) in Marcon et al., 2014).}
}

\note{

}

\value{
  A matrix.
}

\references{
  Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. (2014)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials.
  \emph{Submitted}.
}

\seealso{\code{\link{beed}}, \code{\link{bp}}, \link{beed.confband}}.}

\examples{
v <- as.matrix(0:8, ncol = 1)
Convexity(v, d = 2)
}
\keyword{Convexity, Bernstein polynomials}





