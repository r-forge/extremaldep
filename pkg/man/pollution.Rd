\name{pollution}
\alias{pollution}
\alias{NSN}
\alias{PNN}
\alias{PNS}
\alias{PNNS}
\alias{winterdat}
\alias{Leeds.frechet}

\title{Air quality datasets containing daily maxima of air pollutants (PM10, NO, NO2, 03 and S02)
	recorded in Leeds (U.K.), during five winter seasons (November-Februrary) between 1994 and 1998.}
\description{Contains \eqn{6}{6} datasets: \code{PNS}, \code{PNN}, \code{NSN}, \code{PNNS},
	\code{winterdat} and \code{Leeds.frechet}. }
\details{
	The dataset \code{winterdat} contains \eqn{590}{590} (transformed) observations for 
	each of the five pollutants. Contains \code{NA}s. 
	Outliers have been removed according to Heffernan and Tawn (2004).
	The following datasets have been obtained by applying transformations to \code{winterdat}.
	
	\code{Leeds.frechet} contains \eqn{590}{590} observations corresponding to the
	daily maxima of five air pollutants transformed to unit Frechet scale.
	
	\code{NSN} contains \eqn{100}{100} observations in the \eqn{3}{3}-dimensional
	unit simplex for the daily maxima of nitrogen dioxide (NO2), sulfur dioxide (SO2) 
	and nitrogen oxide (NO).
	
	\code{PNN} contains \eqn{100}{100} observations in the \eqn{3}{3}-dimensional
	unit simplex for the daily maxima of particulate matter (PM10), nitrogen oxide (NO) 
	and nitrogen dioxide (NO2).
	
	\code{PNS} contains \eqn{100}{100} observations in the \eqn{3}{3}-dimensional
	unit simplex for the daily maxima of particulate matter (PM10), nitrogen oxide (NO) 
	and sulfur dioxide (SO2).
	
	\code{PNNS} contains \eqn{100}{100} observations in the \eqn{4}{4}-dimensional
	unit simplex for the daily maxima of particulate matter (PM10), nitrogen oxide (NO), 
	nitrogen dioxide (NO2) and sulfur dioxide (S02).
	
	
	The transformation to unit Frechet margins of the raw data has been considered by
	Cooley et al (2010). Only the \eqn{100}{100} data points with the largest radial components
	were kept.
	}
\source{\url{https://uk-air.defra.gov.uk/}}
\references{

	Cooley, D.,Davis, R. A., and Naveau, P. (2010).
	The pairwise beta distribution: a flexible parametric multivariate model for extremes.
	\emph{Journal of Multivariate Analysis}, \bold{101}, 2103--2117.	
	
	Heffernan, J. E., and Tawn, J. A. (2004).
	A conditional approach for multivariate extreme values.
	\emph{Journal of the Royal Statistical Society, Series B, Methodology}, \bold{66}, 497--546
	
}

\keyword{datasets}