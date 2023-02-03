\name{desn}
\alias{desn}
\alias{pesn}

\title{Univariate extended skew-normal distribution}
\description{Density function, distribution function for the univariate extended skew-normal (\acronym{ESN}) distribution.}

\usage{
desn(x, location=0, scale=1, shape=0, extended=0)
pesn(x, location=0, scale=1, shape=0, extended=0)
}

\arguments{
	\item{x}{quantile. }
	\item{location}{location parameter. \code{0} is the default.}
	\item{scale}{scale parameter; must be positive. \code{1} is the default.}
	\item{shape}{skewness parameter. \code{0} is the default.}
	\item{extended}{extension parameter. \code{0} is the default.}
}

\value{
	density (\code{desn}), probability (\code{pesn}) from the extended skew-normal distribution with given
	\code{location}, \code{scale}, \code{shape} and \code{extended} parameters or from the skew-normal if \code{extended=0}.
	If \code{shape=0} and \code{extended=0} then the normal distribution is recovered.
 }

\references{
Azzalini, A. (1985).
  A class of distributions which includes the normal ones.
  \emph{Scand. J. Statist.} \bold{12}, 171-178.

Azzalini, A. with the collaboration of Capitanio, A. (2014).
 \emph{The Skew-Normal and Related Families}.
 Cambridge University Press, IMS Monographs series.
}

\examples{

dens1 <- desn(x=1, shape=3, extended=2)
dens2 <- desn(x=1, shape=3)
dens3 <- desn(x=1)
dens4 <- dnorm(x=1)
prob1 <- pesn(x=1, shape=3, extended=2)
prob2 <- pesn(x=1, shape=3)
prob3 <- pesn(x=1)
prob4 <- pnorm(q=1)

}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\keyword{distribution}
