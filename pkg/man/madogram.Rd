\name{madogram}
\alias{madogram}
\title{Estimation of the Pickands Dependence Function}
\description{
  Compute a non-parametric estimate of uni- and multi-variate Pickands dependence function, \eqn{A(w)}, based on the madogram.
}
\usage{
  madogram(w, data)
}
\arguments{
  \item{w}{A \eqn{m \times d}{m x d} design matrix (see \bold{Details}).}
  \item{data}{A matrix or data frame with \code{d} columns. \code{d} is the size of the multivariate problem.}
}

\details{
The estimation procedure is based on the madogram as proposed in Marcon et al. (2014). Define the madogram as

\eqn{
\nu(\bold{w}) =
\mathbb{E} \left(\
\bigvee_{i=1,\dots,d}\left \lbrace F^{1/w_i}_{i}\left(X_{i}\right)  \right\rbrace -
\frac{1}{d}\sum_{i=1,\dots,d}F^{1/w_i}_{i}\left(X_{i}\right).
\right),
}{\nu(w) = E[max_{i=1,...,d}[F^{1/w_i}_{i}(X_{i})] - 1/d \sum_{i=1,...,d}F^{1/w_i}_{i}(X_{i}) ],}

with  \eqn{0<w_i<1} and \eqn{w_d= 1 - (w_1 + ... + w_{d-1})}.

For the \eqn{m \times d}{m x d} design matrix \code{w}, each row is a point in the unit
\code{d}-dimensional simplex. 

If \eqn{X} is a \code{d}-dimensional max-stable distributed random vector with exponent measure function \eqn{V} and Pickands dependence function \eqn{A}, then 

\eqn{
\nu(\bold{w})=V(1/w_1,\ldots,1/w_d)/[1+V(1/w_1,\ldots,1/w_d)]-c(\bold{w}),
}
{
\nu(w)=V(1/w_1,\ldots,1/w_d)/[1+V(1/w_1,\ldots,1/w_d)}]-c(w),
}

where \eqn{c(\bold{w})=d^{-1}\sum_{i=1}^{d}{w_i/(1+w_i)}}{c(w)=d^{-1}\sum_{i=1}^{d}{w_i/(1+w_i)}}.
It follows that 

\eqn{
V(1/w_1,\ldots,1/w_d)=\frac{\nu(\bold{w})+c(\bold{w})}{1-\nu(\bold{w})-c(\bold{w})},
}{
V(1/w_1,\ldots,1/w_d)=\frac{\nu(w)+c(w)}{1-\nu(w)-c(\bold{w})},
}
and
\eqn{
A(\bold{w})=\frac{\nu(\bold{w})+c(\bold{w})}{1-\nu(\bold{w})-c(\bold{w})}.
}{
A(w)=\frac{\nu(w)+c(w)}{1-\nu(w)-c(w)}.
}

}

\note{
The madogram-based estimator may not satisfy all constraints required from a Pickands dependence function.
}

\value{
  A numeric vector of estimates.
}

\seealso{\code{\link{beed}}, \code{\link{beed.confband}}}

\examples{
w <- seq(0, 1, length = 100)

data <- rbvevd(50, dep = 0.4, model = "log", mar1 = c(1,1,1))

Amd <- madogram2(w = w, data = data)
Amd.bp <- beed(data = data,x = cbind(w,1-w), d = 2, k = 20, est = 'md', margin ="frechet-exp")
      
plot(w, w, type = "n", xlab = "w1", ylab = "w2", ylim = c(.5,1))
polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = 1, lwd = 1, border = 'grey')
lines(w, Amd, lty = 1, col = 1)
lines(w, Amd.bp$A, lty = 1, col = 2)
}
\keyword{madogram}







