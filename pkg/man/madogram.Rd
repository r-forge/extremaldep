\name{madogram}
\alias{madogram}
\title{Estimation of the Pickands Dependence Function}
\description{
  Compute a non-parametric estimate of the uni- or multi-variate
  Pickands dependence function, \eqn{A(w)}, based on the madogram estimator.
}
\usage{
  madogram(w, data)
}
\arguments{
  \item{w}{A \eqn{(m \times d)}{(m x d)} design matrix (see \bold{Details}).}
  \item{data}{A \eqn{(n \times d)}{(n x d)} matrix of data or data frame with \code{d} columns. \code{d} is
    the numer of variables and \code{n} is the number of replications.}
}

\details{
The estimation procedure is based on the madogram as proposed in Marcon et al. (2014). The madogram is defined by

\eqn{
\nu(\bold{w}) =
\mathbb{E} \left(\
\bigvee_{i=1,\dots,d}\left \lbrace F^{1/w_i}_{i}\left(X_{i}\right)  \right\rbrace -
\frac{1}{d}\sum_{i=1,\dots,d}F^{1/w_i}_{i}\left(X_{i}\right).
\right),
}{\nu(w) = E[max_{i=1,...,d}(F^{1/w_i}_{i}(X_{i})) - 1/d \sum_{i=1,...,d}F^{1/w_i}_{i}(X_{i})],}

where \eqn{0<w_i<1}{0<w_i<1} and \eqn{w_d=1-(w_1+\ldots+w_{d-1})}{w_d=1-(w_1+...+w_{d-1})}.

Each row of the design matrix \code{w} is a point in the unit
\code{d}-dimensional simplex. 

If \eqn{X}{X} is a \code{d}-dimensional max-stable distributed random vector, with exponent measure function \eqn{V(\bold{x})}{V(x)} and Pickands dependence function \eqn{A(\bold{w})}{A(w)}, then 

\eqn{\nu(\bold{w}=V(1/w_1,\ldots,1/w_d)/(1+V(1/w_1,\ldots,1/w_d))-c(\bold{w}),}{\nu(w) = V(1/w_1,...,1/w_d)/(1+V(1/w_1,...,1/w_d))-c(w),}
where \eqn{c(\bold{w})=d^{-1}\sum_{i=1}^{d}{w_i/(1+w_i)}}{c(w)=d^{-1}\sum_{i=1}^{d}{w_i/(1+w_i)}}.

From this, it follows that 

\eqn{
V(1/w_1,\ldots,1/w_d)=\frac{\nu(\bold{w})+c(\bold{w})}{1-\nu(\bold{w})-c(\bold{w})},
}{
V(1/w_1,...,1/w_d)=(\nu(w)+c(w))/(1-\nu(w)-c(w)),
}

and

\eqn{
A(\bold{w})=\frac{\nu(\bold{w})+c(\bold{w})}{1-\nu(\bold{w})-c(\bold{w})}.
}{
A(w)=(\nu(w)+c(w))/(1-\nu(w)-c(w)).
}
}

\note{
The madogram-based estimator does not satisfy the required constraints in order to provide a valid Pickands dependence function.
}

\value{
  A numeric vector of estimates.
}

\references{
  Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. (2014)
  Multivariate Nonparametric Estimation of the Pickands Dependence
  Function using Bernstein Polynomials,
  \emph{Submitted}.
  
  Naveau, P., Guillou, A., Cooley, D., Diebolt, J. (2009)
  Modelling pairwise dependence of maxima in space,
  Biometrika, \bold{96(1)}, 1-17.
}

\seealso{\code{\link{beed}}, \code{\link{beed.confband}}}

\examples{
w <- seq(0, 1, length = 100)

data <- rbvevd(50, dep = 0.4, model = 'log', mar1 = c(1,1,1))

Amd <- madogram(cbind(w,1-w), data)
Amd.bp <- beed(data, cbind(w,1-w), 2, 'md', 'emp', 20, plot=TRUE)
      
lines(w, Amd, lty = 1, col = 2)
}
\keyword{madogram}
