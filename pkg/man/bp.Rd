\name{bp}
\alias{bp}
\title{Univariate and Multivariate Bernstein Polynomials}
\description{
  Compute the uni- or multi-variate Bernstein polynomials of oreder \code{k}.
}
\usage{
  bp(x, v, k)
}
\arguments{
  \item{x}{A \eqn{m \times d}{m x d} design matrix (see \bold{Details}).}
  \item{v}{A matrix of the index set of \eqn{\beta}{beta} coefficients.}
  \item{k}{An integer parameter for the order of the Bernstein polynomials.}
}

\details{
For the \eqn{m \times d}{m x d} design matrix \code{x}, each row is a point in the unit \code{d}-dimensional simplex,

\eqn{
S_d := \left\{ (w_1,\ldots, w_{d-1}) \in [0,1]^{d-1}: \sum_{i=1}^{d-1} w_i \leq 1 \right\}.
}{
S_d := { (w_1,..., w_{d-1}) in [0,1]^{d-1}: \sum_{i=1}^{d-1} w_i <= 1}.
}

We define the Bernstein basis polynomial of degree \code{k} as the following continuous function with values in \eqn{[0,1]},

\eqn{
b(\bold{w};k):= \binom{k}{\bold{\alpha}_l} \;\bold{w}^{\bold{\alpha}_l} \;(1-\Sigma_{\bold{w}})^{k-\Sigma_{\bold{\alpha}_l}},
}{
b(w;k):= \binom{k}{\alpha_l} w^{\alpha_l} (1-\Sigma_{w})^{k-\Sigma_{\alpha_l}},
}

where \eqn{\bold{w}}{w} are points of the \code{d}-dimensional unit simplex \eqn{S_d} and, 

\eqn{
\binom{k}{\bold{\alpha}_l} = \frac{k!}{\bold{\alpha}_l! (k-\Sigma_{\bold{\alpha}_l})!}, \qquad \bold{\alpha}_l! = j_1! \cdots j_{d-1}\,!,  \qquad \bold{w}^{\bold{\alpha}_l} = \prod_{i=1}^{d-1} w_i^{j_i}.
}{
\binom{k}{\alpha_l} = \frac{k!}{\alpha_l! (k-\Sigma_{\alpha_l})!}, 
\alpha_l! = j_1!...j_{d-1}!, w^{\alpha_l} = \prod_{i=1}^{d-1} w_i^{j_i}.
}

Therefore, the Bernstein-B\'{e}zier polynomial representation of a function \eqn{f} is given by

\eqn{
B_f(\bold{w};k) = \sum_{l \in L_k} \beta_l b_l (\bold{w};k).
}{
B_f(w;k) = \sum_{l \in L_k} beta_l b_l (w;k).
}
}

%\note{
  
  \value{
    A vector/matrix of Bernstein polynomial basis.
  }
  
  \references{
Marcon, G., Padoan, S.A., Naveau, P., Muliere, P. (2014)
Multivariate Nonparametric Estimation of the Pickands Dependence
Function using Bernstein Polynomials.
\emph{Submitted}.
    
Petrone, S. (2004) 
On the role of mixtures in bayesian nonparametrics. 
\emph{Atti della XLII Riunione Scientifica della Societa Italiana di Statistica 1}, 257--268.

Li, F. J. (2011) 
Interpolation and convergence of Bernstein-B\'{e}zier coefficients.
\emph{Acta Mathematica Sinica, English Series} \bold{27 (9)}, 1769--1782.
  }
  
\seealso{\code{\link{beed}}, \code{\link{beed.confband}}.}
  
\examples{
w <- seq(0,1, length=100)    
data <- rbvevd(50, dep = 0.4, model = "log", mar1 = c(1,1,1))
    
k = 20
beta.index <- expand.grid(0:k, 0:k)
beta.index.sum <- apply(beta.index, 1, sum)
restr <- which(beta.index.sum <= k)
v <- as.matrix(beta.index[restr, ncol(beta.index):1])
 
bpbasis <- bp(x = cbind(w, 1-w), v = v, k = k)
Amd <- beed(data = data, x = cbind(w, 1-w), d = 2, k = 20, 
            est = 'md', margin ="frechet-exp")
}

\keyword{bernstein polynomials}
  
  
  
  
  
  