\name{ellipse}
\alias{ellipse}

\title{Level sets for bivariate normal, student-t and skew-normal distributions probability densities.}
\description{Level sets of the bivariate normal, student-t and skew-normal distributions probability densities for a given probability.}

\usage{
	ellipse(center=c(0,0), alpha=c(0,0), sigma=diag(2), df=1,
	prob=0.01, npoints=250, pos=FALSE)
}

\arguments{
	\item{center}{A vector of length 2 corresponding to the location of the distribution.}
	\item{alpha}{A vector of length 2 corresponding to the skewness of the skew-normal distribution.}
	\item{sigma}{A 2 by 2 variance-covariance matrix.}
	\item{df}{An integer corresponding to the degree of freedom of the student-t distribution.}
	\item{prob}{The probability level. See \code{details}}
	\item{npoints}{The maximum number of points at which it is evaluated.}
	\item{pos}{If \code{pos=TRUE} only the region on the positive quadrant is kept.}
}

\value{
	Returns a bivariate vector of \eqn{250}{250} rows if \code{pos=FALSE}, and half otherwise.
}

\details{
	The Level sets are defined as

	\eqn{R(f_\alpha)=\{ x: f(x) \geq f_\alpha \} }{R(f_\alpha)=\{ x: f(x) \geq f_\alpha \}}

	where \eqn{f_\alpha}{f_\alpha} is the largest constant such that

	\eqn{P(X \in R(f_\alpha)) \geq 1-\alpha}{P(X \in R(f_\alpha)) \geq 1-\alpha}.
	Here we consider \eqn{f(x)}{f(x)} to be the bivariate normal, student-t or skew-normal density.

}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\examples{

library(mvtnorm)

# Data simulation (Bivariate-t on positive quadrant)
rho <- 0.5
sigma <- matrix(c(1,rho,rho,1), ncol=2)
df <- 2

set.seed(101)
n <- 1500
data <- rmvt(5*n, sigma=sigma, df=df)
data <- data[data[,1]>0 & data[,2]>0, ]
data <- data[1:n, ]

P <- c(1/750, 1/1500, 1/3000)

ell1 <- ellipse(prob=1-P[1], sigma=sigma, df=df, pos=TRUE)
ell2 <- ellipse(prob=1-P[2], sigma=sigma, df=df, pos=TRUE)
ell3 <- ellipse(prob=1-P[3], sigma=sigma, df=df, pos=TRUE)

plot(data, xlim=c(0,max(data[,1],ell1[,1],ell2[,1],ell3[,1])),
     ylim=c(0,max(data[,2],ell1[,2],ell2[,2],ell3[,2])), pch=19)
points(ell1, type="l", lwd=2, lty=1)
points(ell2, type="l", lwd=2, lty=1)
points(ell3, type="l", lwd=2, lty=1)

}

\keyword{models}
