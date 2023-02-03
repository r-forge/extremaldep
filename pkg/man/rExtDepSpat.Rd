\name{rExtDepSpat}
\alias{rExtDepSpat}

\title{Random generation of max-stable processes}
\description{This function generates realisations from a max-stable process.}

\usage{
rExtDepSpat(n, coord, model="SCH", cov.mod = "whitmat", grid = FALSE, 
            control = list(), cholsky = TRUE, ...)
}

\arguments{
	\item{n}{An integer indictaing the number of observations.}
	\item{coord}{A vector or matrix corresponding to the coordinates of locations where the processes is simulated. Each row corresponds to a location.}
	\item{model}{A character string indicating the max-stable model. See \code{details}.}
	\item{cov.mod}{A character string indicating the correlation function function. See \code{details}.}
	\item{grid}{A logical value; \code{TRUE} for coordinates on a grid.}
	\item{control}{A named list with arguments \code{nlines} giving the number of lines of the TBM simulation, \code{method} a character string specifying the name of the simulation method and \code{uBound} the uniform upper bound. Note that \code{method} must take value \code{'exact'}, \code{'tbm'} or \code{'circ'}. See \code{details}.}
	\item{cholsky}{A logical value; if \code{TRUE} a Cholesky decomposition is performed. Used for the extremal-t and extremal skew-t models when \code{control=list(method='exact')}. }
	\item{...}{The parameters of the max-stable model. See \code{details}.}
}

\details{
  This function extends the \code{rmaxstab} function from the \code{SpatialExtremes} package in two ways:
  \itemize{
  \item{1. }{The \bold{extremal skew-t} model is included.}
  \item{2. }{The function returns the hitting scenarios, i.e. the index of which 'storm' (or process) led to the maximum value for each location and observation.}
  }

  The max-stable models available in this procedure and the specifics  are:
  \itemize{
    \item{\bold{Smith model}: }{when \code{model='SMI'}, does not require \code{cov.mod}. If \code{coord} is univariate then \code{var} needs to be specified and for higher dimensions covariance parameters should be provided such as \code{cov11}, \code{cov12}, \code{cov22}, etc. }
    \item{\bold{Schlather model}: }{when \code{model='SCH'}, requires \code{cov.mod='whitmat'}, \code{'cauchy'}, \code{'powexp'} or \code{'bessel'} depending on the correlation family. Parameters \code{'nugget'}, \code{'range'} and \code{'smooth'} should be specified.}
    \item{\bold{Extremal-t model}: }{when \code{model='ET'}, requires \code{cov.mod='whitmat'}, \code{'cauchy'}, \code{'powexp'} or \code{'bessel'} depending on the correlation family. Parameters \code{'nugget'}, \code{'range'}, \code{'smooth'} and \code{'DoF'} should be specified.}
    \item{\bold{Extremal skew-t model}: }{when \code{model='EST'}, requires \code{cov.mod='whitmat'}, \code{'cauchy'}, \code{'powexp'} or \code{'bessel'} depending on the correlation family. Parameters \code{'nugget'}, \code{'range'}, \code{'smooth'}, \code{'DoF'}, \code{'alpha'} (a vector of length \eqn{3}{3}) and \code{'acov1'} and \code{'acov2'} (both vector of length the number of locations) should be specified. The skewness vector is defined as \eqn{\alpha = \alpha_0 + \alpha_1 \textrm{acov1} + \alpha_2 \textrm{acov2} }{\alpha = \alpha_0 + \alpha_1 \textrm{acov1} + \alpha_2 \textrm{acov2}}.}
    \item{\bold{Geometric gaussian model}: }{when \code{model='GG'}, requires \code{cov.mod='whitmat'}, \code{'cauchy'}, \code{'powexp'} or \code{'bessel'} depending on the correlation family. Parameters \code{'sig2'}, \code{'nugget'}, \code{'range'} and \code{'smooth'} should be specified.}
    \item{\bold{Brown-Resnick model}: }{when \code{model='BR'}, does not require \code{cov.mod}. Parameters \code{'range'} and \code{'smooth'} should be specified.}  
  }
  
  For the argument \code{control}, details of the list components are as follows:
  \itemize{
    \item{\code{method}}{is \code{NULL} by default, meaning that the function tries to find the most appropriate simulation technique. Current simulation techniques are a direct approach, i.e. Cholesky decomposition of the covariance matrix, the turning bands and the circular embedding methods. Note that for the extremal skew-t model it can only take value \code{'exact'} or \code{'direct'};}
    \item{\code{nlines}}{ if \code{NULL} then it is set to \eqn{1000}{1000};}
    \item{\code{uBound}}{ if \code{NULL} then it is set to reasonable values - for example \eqn{3.5}{3.5} for the Schlather model.}
  }
}

\value{
  A list made of
  \itemize{
    \item{\code{vals}: }{A \eqn{(n \times d)}{(n x d)} matrix containing \eqn{n}{n} observations at \eqn{d}{d} locations, from the specified max-stable model.}
    \item{\code{hits}: }{A \eqn{(n \times d)}{(n x d)} matrix containing the hitting scenarios for each observations. On each row, elements with the same integer value indicate that the maxima at these two locations is coming from the same 'storm' or process.}
  }
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\references{
	Beranger, B., Stephenson, A. G. and Sisson, S.A. (2021)
  High-dimensional inference using the extremal skew-t process
  \emph{Extremes}, \bold{24}, 653-685.
}

\seealso{\code{\link{fExtDepSpat}} }

\examples{
# Generate some locations
set.seed(1)
lat <- lon <- seq(from=-5, to=5, length=20)
sites <- as.matrix(expand.grid(lat,lon))

# Example using the extremal-t
set.seed(2)
z <- rExtDepSpat(1, sites, model="ET", cov.mod="powexp", DoF=1, 
                 nugget=0, range=3, smooth=1.5, 
                 control=list(method="exact"))
fields::image.plot(lat, lon, matrix(z$vals,ncol=20) )

# Example using the extremal skew-t
if(interactive()){
set.seed(3)
z2 <- rExtDepSpat(1, sites, model="EST", cov.mod="powexp", DoF=5, 
                  nugget=0, range=3, smooth=1.5, alpha=c(0,5,5), 
                  acov1=sites[,1], acov2=sites[,2], 
                  control=list(method="exact"))
fields::image.plot(lat, lon, matrix(z2$vals,ncol=20) )
}
}

\keyword{models}
