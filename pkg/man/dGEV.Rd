\name{dGEV}
\alias{dGEV}
\alias{pGEV}
\alias{qGEV}

\title{The Generalized Extreme Value Distribution}
\description{Density, distribution and quantile function for the Generalized Extreme Value (\acronym{GEV}) distribution.}

\usage{
dGEV(x, loc, scale, shape, log=FALSE)
pGEV(q, loc, scale, shape, lower.tail=TRUE)
qGEV(p, loc, scale, shape)
}

\arguments{
	\item{x, q}{vector of quantiles. }
	\item{p}{vector of probabilities.}
	\item{loc}{vector of locations.}
	\item{scale}{vector of scales.}
	\item{shape}{vector of shapes.}
	\item{log}{A logical value; if \code{TRUE} returns the log density.}
	\item{lower.tail}{A logical value; if \code{TRUE} probabilities are \eqn{P\left(X \leq x \right)}{P(X\leq x)}, otherwise \eqn{P\left(X > x \right)}{P(X > x)}.}
}

\details{
  The GEV distribution has density
  \deqn{f(x; \mu, \sigma, \xi) = \exp \left\{ -\left[ 1 + \xi \left( \frac{x-\mu}{\sigma} \right)\right]_+^{-1/\xi}\right\}}{f(x; \mu,\sigma,\xi) = exp \{ -[1 + \xi ( (x-\mu)/\sigma )]_+^(-1/\xi) \} }
}

\value{
	density (\code{dGEV}), distribution function (\code{pGEV}) and quantile function (\code{qGEV}) from the Generalized Extreme Value distirbution with given
	\code{location}, \code{scale} and \code{shape}.
}

\seealso{ \code{\link{fGEV}} }

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\examples{
# Densities
dGEV(x=1, loc=1, scale=1, shape=1)
dGEV(x=c(0.2, 0.5), loc=1, scale=1, shape=c(0,0.3))

# Probabilities
pGEV(q=1, loc=1, scale=1, shape=1, lower.tail=FALSE)
pGEV(q=c(0.2, 0.5), loc=1, scale=1, shape=c(0,0.3))

# Quantiles
qGEV(p=0.5, loc=1, scale=1, shape=1)
qGEV(p=c(0.2, 0.5), loc=1, scale=1, shape=c(0,0.3))
}



\keyword{distribution}
