\name{heat}
\alias{heat}
\alias{geogmel}
\alias{heatdata}
\alias{locgrid}
\alias{mellat}
\alias{mellon}
\alias{scalegrid}
\alias{shapegrid}

\title{Summer temperature maxima in Melbourne, Australia between 1961 and 2010.}

\description{The dataset corresponds to the summer maxima taken over the period from August to April inclusive, recorded between 1961 and 2010 at 90 stations on a 0.15 degree grid in a 9 by 10 formation.}

\details{The first maximum is taken over the August 1961 to April 1962 period, and the last maximum is taken over the August 2010 to April 2011 period. The object \code{heatdata} contains the core of the data:
  \itemize{
    \item{\code{vals}: }{A \eqn{50*90}{50*90} matrix containing the \eqn{50}{50} summer maxima at the \eqn{90}{90} locations.}
    \item{\code{sitesLL: }A \eqn{90*2}{90*2} matrix containing the sites locations in Latitude-Longitude, recentred (means have been substracted).}
    \item{\code{sitesEN: }A \eqn{90*2}{90*2} matrix containing the sites locations in Eastings-Northings, recentred (means have been substracted).}
    \item{\code{vals}: }{A \eqn{50*90}{50*90} matrix containing integers indicating the ``heatwave'' number of each of the \eqn{50}{50} summer maxima at all \eqn{90}{90} locations. Locations on the same row with the same integer indicates that they were obtained from the same heatwave. Heatwaves are defined over a three day window.}
    \item{\code{sitesLLO: }A \eqn{90*2}{90*2} matrix containing the sites locations in Latitude-Longitude, on the original scale.}
    \item{\code{sitesENO: }A \eqn{90*2}{90*2} matrix containing the sites locations in Eastings-Northings, on the original scale.}
    \item{\code{ufvals}: }{A \eqn{50*90}{50*90} matrix containing the \eqn{50}{50} summer maxima at the \eqn{90}{90} locations, on the unit Frechet scale.}
  }
  
  Standardisation to unit Frechet is performed as in Beranger et al. (2021) by fitting the GEV distribution marginally using unconstrained location and shape parameters and the shape parameter to be a linear function of eastings and northings in 100 kilometre units. The resulting estimates are given in the objects \code{locgrid}, \code{scalegrid} and \code{shapegrid}, which are \eqn{10*9}{10*9} matrices.
  
  Details about the study region are given in \code{geogmel}, \code{mellat} and \code{mellon}. \code{geogmel} is a large SpatialPolygonsDataFrame of class s4 that was obtained using the \code{sp} package. \code{mellat} and \code{mellon} are vectors of length \eqn{10}{10} and \eqn{11}{11} which give the latitude and longitude coordinates of the grid.
}

\references{
	Beranger, B., Stephenson, A. G. and Sisson, S.A. (2021)
  High-dimensional inference using the extremal skew-t process
  \emph{Extremes}, \bold{24}, 653-685.
}

\examples{
###############################
### Plot the study region   ###
### using the rgeos library ###
###############################

if(interactive()){  

plot(rgeos::gSimplify(geogmel, 0.01), lwd=1, col=grey(0.75))
points(heatdata$sitesLLO, pch=16)

}
}

\keyword{datasets}
