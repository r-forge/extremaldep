\name{PNNS}
\alias{PNNS}
\alias{print.PNNS}

\title{A \eqn{4}{4}-dimensional dataset with unit Frechet margins.}
\description{The data contains \eqn{100}{100} observations in the \eqn{4}{4}-dimensional
	unit simplex. These transformed observations correspond to the daily maxima of 
	particulate matter (PM10), nitrogen oxide (NO), nitrogen dioxide (NO2) and sulfur dioxide (S02)
	recorded in the city of Leeds during the winter season (November-February) between 1994 and 1998.}
\details{The transformation to unit Frechet margins of the raw data has been considered by
	Cooley et al (2010). Only the \eqn{100}{100} data points with the largest radial components
	were kept.}	
\format{A \eqn{100*4}{100*4} matrix.}
\source{\url{http://airquality.co.uk}}
\references{
	Cooley, D.,Davis, R. A., and Naveau, P. (2010).
	The pairwise beta distribution: a flexible parametric multivariate model for extremes.
	\emph{Journal of Multivariate Analysis}, \bold{101}, 2103--2117.	
}
\seealso{
	\link{PNS},
	\link{NSN},
	\link{PNN}.
}		