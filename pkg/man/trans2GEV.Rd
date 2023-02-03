\name{trans2GEV}
\alias{trans2GEV}

\title{Transformation to GEV distribution}
\description{Transformation of marginal distribution from unit Frechet to GEV}

\usage{
trans2GEV(data, pars)
}

\arguments{
	\item{data}{A vector of length \eqn{n}{n} or a \eqn{(n \times p)}{(n x p)} matrix representing the data on its original scale.}
	\item{pars}{A \eqn{(1 \times 3)}{(1 x 3)} vector or a \eqn{(p \times 3)}{(p x 3)} matrix of marginal GEV parameters.}
}

\details{
 The transformation function is \eqn{\left(x^{\xi} -1 \right) \frac{\sigma}{\xi}+\mu}{(x^\xi -1)\sigma/\xi + \mu} if \eqn{\xi\neq 0}{\xi \neq 0}, and \eqn{x^{-1}\sigma+\mu}{\sigma/x + \mu} if \eqn{\xi=0}{\xi=0}. 
}

\value{
	An object of the same format and dimensions as \code{data}.   
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\seealso{\code{\link{trans2UFrechet}} }

\examples{
data(pollution)
pars <- fGEV(Leeds.frechet[,1])$est

par_new <- c(2, 1.5, 0.5)
data_new <- trans2GEV(Leeds.frechet[,1], pars=par_new)

fGEV(data_new)
}

\keyword{models}
