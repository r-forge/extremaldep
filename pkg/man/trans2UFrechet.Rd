\name{trans2UFrechet}
\alias{trans2UFrechet}

\title{Transformation to unit Frechet distribution}
\description{Empirical and parametric transformation of a dataset to unit Frechet marginal distribution}

\usage{
trans2UFrechet(data, pars, type="Empirical")
}

\arguments{
	\item{data}{A vector of length \eqn{n}{n} or a \eqn{(n \times p)}{(n x p)} matrix representing the data on its original scale.}
	\item{pars}{A \eqn{(1 \times 3)}{(1 x 3)} vector or a \eqn{(p \times 3)}{(p x 3)} matrix of marginal GEV parameters. Required when \code{type="GEV"}.}
	\item{type}{A character string indicating the type of transformation. Can take value \code{"Empirical"} or \code{"GEV"}.}
}

\details{
  When \code{type="Empirical"}, the transformation function is \eqn{t(x)=-1/\log(F_{\textrm{emp}}(x))}{t(x = -1/log(F_{emp}(x)))} where \eqn{F_{\textrm{emp}}(x)}{F_{emp}(x)} denotes the empirical cumulative distribution.
  
  When \code{type="GEV"}, the transformation function is \eqn{\left(1+\xi \frac{x-\mu}{\sigma}\right)^{1/\xi}}{(1+\xi \frac{x-\mu}{\sigma})^{1/\xi}} if \eqn{\xi \neq 0}{\xi \neq 0}, \eqn{\left( \frac{x-\mu}{\sigma}\right)^{-1}}{\sigma / (x-\mu)} if \eqn{\xi=0}{\xi=0}. If the argument \code{pars} is missing then a GEV is fitted on the columns of \code{data} using the \code{\link{fGEV}} function.
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

\seealso{\code{\link{trans2GEV}}, \code{\link{fGEV}} }

\examples{
data(MilanPollution)

pars <- fGEV(Milan.winter$PM10)$est
pars

data_uf <- trans2UFrechet(data=Milan.winter$PM10, pars=pars, 
                          type="GEV")
fGEV(data_uf)$est

data_uf2 <- trans2UFrechet(data=Milan.winter$PM10, 
                           type="Empirical")
fGEV(data_uf2)$est
}

\keyword{models}
