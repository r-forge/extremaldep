\name{Leeds.frechet}
\alias{Leeds.frechet}
\alias{print.Leeds.frechet}

\title{\eqn{5}{5}-dimensional dataset with unit Frechet margins.}

\description{The dataset contains \eqn{590}{590} (transformed) daily maxima of five 
	air pollutants recorded in Leeds (U.K.) during five winter seasons (November-Februrary)
	between 1994 and 1998. Contains \code{NA}'s. Marginal transformation to unit Frechet was perfomed
	by Cooley et al. (2010). Outliers have been removed according to Heffernan and Tawn (2004).}
  
\details{The transformation to unit Frechet margins of the raw data has been considered by
	Cooley et al (2010). Only the \eqn{100}{100} data points with the largest radial components
	were kept.}
  
\format{A \eqn{590*5}{590*5} matrix.}

\source{\url{http://airquality.co.uk}}

\references{
	Cooley, D.,Davis, R. A., and Naveau, P. (2010).
	The pairwise beta distribution: a flexible parametric multivariate model for extremes.
	\emph{Journal of Multivariate Analysis}, \bold{101}, 2103--2117.	
	
	Heffernan, J. E., and Tawn, J. A. (2004).
	A conditional approach for multivariate extreme values.
	\emph{Journal of the Royal Statistical Society, Series B, Methodology}, \bold{66}, 497--546
	
}
