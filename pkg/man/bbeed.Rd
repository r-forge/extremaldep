\name{bbeed}
\alias{bbeed}
\title{Bayesian Estimation of Extremal Dependence}
\description{
Nonparametric Bayesian estimation of the extremal dependence function (angular measure and Pickands dependence function) for two-dimensional data on the basis of Bernstein polynomials representations.
}
\usage{
bbeed(data, pm0, param0, k0, hyperparam, nsim, nk = 70,
      prior.k = c('pois','nbinom'), prior.pm = c('unif','beta', 'null'),
      lik = TRUE)
}
\arguments{
  \item{data}{\eqn{(n \times d)}{(n x 2)} matrix of component-wise maxima (\code{n} is the number of replications).}
  \item{pm0}{list; the initial state of the Markov chain of point masses \eqn{p_0}{p0} and \eqn{p_1}{p1}. Randomly generated if missing.}
  \item{param0}{real vector in (0,1) of length \code{k0}+1; the initial state of the Markov chain of \eqn{\eta}{\eta} coefficients (see \bold{Details}). Randomly generated if missing.}
  \item{k0}{postive integer indicating the initial state of the Markov chain of the Bernstein polynomial's degree. Randomly generated if missing.}
  \item{hyperparam}{list containing the hyperparameters (see \bold{Details}).}
  \item{nsim}{postive integer indicating the number of the iterations of the chain.}
  \item{nk}{positive integer, maximum value of \eqn{k} that can be reached (\code{nk = 70} by default).}
  \item{prior.k}{string, indicating the prior distribution on the degree of the Bernstein polynomial, \eqn{k}. If missing, \code{prior.k = "pois"} by default.}
  \item{prior.pm}{string, indicating the prior distribution on the point masses \eqn{p_0}{p0} and \eqn{p_1}{p1}. If missing, \code{prior.pm = "unif"} by default.}
  \item{lik}{logical; if \code{FALSE}, an approximation of the likelihood, by means of the angular measure in Bernstein form is used (\code{TRUE} by default). }
}

\details{

The hyperparameters and starting values may not be specified and in this case a warning message will be printed and default values are set. In particular if \code{hyperparam} is missing we have

\code{
	hyperparam <- NULL
  hyperparam$a.unif <- 0
  hyperparam$b.unif <- .5
  hyperparam$a.beta <- c(.8,.8)
  hyperparam$b.beta <- c(5,5)
  mu.pois <- hyperparam$mu.pois <- 4
  mu.nbinom <- hyperparam$mu.nbinom <- 4
  var.nbinom <- 8
  pnb <- hyperparam$pnb <- mu.nbinom/var.nbinom
  rnb <- hyperparam$rnb <- mu.nbinom^2 / (var.nbinom-mu.nbinom)
}

  The routine returns an estimate of the Pickands dependence function using the Bernstein polynomials approximation proposed in Marcon et al. (2016).
  The method is based on a preliminary empirical estimate of the Pickands dependence function.
  If you do not provide such an estimate, this is computed by the routine. In this case, you can select one of the empirical methods
  available. \code{est = "ht"} refers to the Hall-Tajvidi estimator (Hall and Tajvidi 2000).
  With \code{est = "cfg"} the method proposed by Caperaa et al. (1997) is considered. Note that in the multivariate case the adjusted  version of Gudendorf and Segers (2011) is used. Finally, with \code{est = "md"} the estimate is based on the madogram defined in Marcon et al. (2016).

  For each row of the \eqn{(m \times d)}{(m x d)} design matrix \code{x} is a point in the unit \code{d}-dimensional simplex,

\eqn{
S_d := \left\{ (w_1,\ldots, w_d) \in [0,1]^{d}: \sum_{i=1}^{d} w_i = 1 \right\}.
}{
S_d := { (w_1,..., w_d) in [0,1]^{d}: \sum_{i=1}^{d} w_i = 1}.
}

  With this "regularization" method, the final estimate satisfies the neccessary conditions in order to be a Pickands dependence function.

  \eqn{A(\bold{w}) = \sum_{\bold{\alpha} \in \Gamma_k} \beta_{\bold{\alpha}} b_{\bold{\alpha}} (\bold{w};k).}{A(w) = \sum_{\alpha in \Gamma_k} \beta_\alpha b_\alpha (w;k).}

  The estimates are obtained by solving an optimization quadratic problem subject to the constraints. The latter are represented by the following conditions:
  \eqn{A(e_i)=1; \max(w_i)\leq A(w) \leq 1; \forall i=1,\ldots,d;}{A(e_i)=1; \max(w_i)<=A(w)<=1; for all i=1,...,d;} (convexity).

  The order of polynomial \code{k} controls the smoothness of the estimate. The higher \code{k} is, the smoother the final estimate is.
  Higher values are better with strong dependence (e. g. \code{k=23}), whereas small values (e.g. \code{k=6} or \code{k=10}) are enough with mild or weak dependence.

An empirical transformation of the marginals is performed when \code{margin="emp"}. A max-likelihood fitting of the GEV distributions is implemented when \code{margin="est"}. Otherwise it refers to marginal parametric GEV theorethical distributions (\code{margin = "exp", "frechet", "gumbel"}).

If \code{bounds = TRUE}, a modification can be implemented to satisfy
\eqn{\max(w,1-w) \leq A_n(w) \leq 1 \forall 0 \leq w \leq 1}{\max(w,1-w) <= A_n(w) <= 1 for all 0 <= w <= 1}:

\eqn{ A_n(w) = \min(1, \max{A_n(w), w, 1-w}). }{ A_n(w) = min(1, max{A_n(w), w, 1-w}). }
}

\note{
The number of coefficients depends on both the order of polynomial \code{k} and the dimension \code{d}. The number of parameters is \eqn{d-1 \choose k+d-1}{choose(k+d-1,d-1)}.

The size of the vector \code{beta} must be compatible with the polynomial order \code{k} chosen.

With the estimated polynomial coefficients, the extremal coefficient, i.e. \eqn{d*A(1/d,\ldots,1/d)}{d*A(1/d,...,1/d)} is computed.
}

\value{

return(list(pm=spm, eta=seta, k=sk, accepted=accepted/nsim,
              prior = list(hyperparam=hyperparam, k=prior.k, pm=prior.pm)))


  \item{beta}{vector of polynomial coefficients}
  \item{A}{numeric vector of the estimated Pickands dependence function \eqn{A}}
  \item{Anonconvex}{preliminary non-convex function}
  \item{extind}{extremal index}
}

\references{
  Marcon G., Padoan, S.A. and Antoniano-Villalobos I. (2016) Bayesian Inference for the Extremal Dependence. \emph{Electronic Journal of Statistics}, \bold{10}(2), 3310-3337.
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.unibocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{http://www.borisberanger.com};
	Giulia Marcon, \email{giuliamarcongm@gmail.com}
}

\seealso{\code{\link{summary.bbeed}}}

\examples{
\donttest{

# This reproduces some of the results showed in Fig. 1 (Marcon, 2016).
set.seed(1890)
data <- evd::rbvevd(n=100, dep=.6, asy=c(0.8,0.3), model="alog", mar1=c(1,1,1))

nsim = 500000
burn = 400000

mu.nbinom = 3.2
var.nbinom = 4.48
hyperparam <- list(a.unif=0, b.unif=.5, mu.nbinom=mu.nbinom, var.nbinom=var.nbinom)
k0 = 5
pm0 = list(p0=0.06573614, p1=0.3752118)
eta0 = ExtremalDep:::rcoef(k0, pm0)

mcmc <- bbeed(data, pm0, eta0, k0, hyperparam, nsim,
              prior.k = "nbinom", prior.pm = "unif")

w <- seq(0.001, .999, length=100)
summary.mcmc <- summary.bbeed(w, mcmc, burn, nsim, plot=TRUE)
}
}

\keyword{Bayesian Nonparametric}
