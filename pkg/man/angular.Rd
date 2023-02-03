\name{angular}
\alias{angular}

\title{Estimation of the angular density, angular measure and random generation from the angular distribution.}

\description{Empirical estimation to the Pickands dependence function, the angular density, the angular measure and random generation of samples from the estimated angular density.}

\usage{
angular(data, model, n, dep, asy, alpha, beta, df, seed, k, nsim, plot=TRUE, nw=100)
}

\arguments{
	\item{data}{The dataset in vector form}
	\item{model}{The specified model; a character string. Must be either \code{"log"}, \code{"alog"}, \code{"hr"}, \code{"neglog"}, \code{"aneglog"}, \code{"bilog"}, \code{"negbilog"}, \code{"ct"}, \code{"amix"} or \code{"Extremalt"} for the logistic, asymmetric logistic, Husler-Reiss, negative logistic, asymmetric negative logistic, bilogistic, negative bilogistic, Coles-Tawn, asymetric mixed and Extremal-t models respectively.}
	\item{n}{The number of random generations from the \code{model}. Required if \code{data=NULL}.}
	\item{dep}{The dependence parameter for the \code{model}.}
	\item{asy}{A vector of length two, containing the two asymmetry parameters for the asymmetric logistic (\code{model='alog'}) and asymmetric negative logistic models (\code{model='aneglog'}).}
	\item{alpha,beta}{Alpha and beta parameters for the bilogistic, negative logistic, Coles-Tawn and asymmetric mixed models.}
	\item{df}{The degree of freedom for the extremal-t model.}
	\item{seed}{The seed for the data generation. Required if \code{data=NULL}.}
	\item{k}{The polynomial order.}
	\item{nsim}{The number of generations from the estimated angular density.}
	\item{plot}{If \code{TRUE}, the fitted angular density, histogram of the generated observations from the angular density and the true angular density (if \code{model} is specified) are displayed.}
	\item{nw}{The number of points at which the estimated functions are evaluated}
}

\details{
	See Marcon et al. (2017).
}

\value{
	Returns a list which contains \code{model}, \code{n}, \code{dep}, \code{data}, \code{Aest} the estimated pickands dependence function, \code{hest} the estimated angular density, \code{Hest} the estimated angular measure, \code{p0} and \code{p1} the point masses at the edge of the simplex, \code{wsim} the simulated sample from the angular density and \code{Atrue} and \code{htrue} the true Pickand dependence function and angular density (if \code{model} is specified).
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
	Giulia Marcon, \email{giuliamarcongm@gmail.com}
}

\references{
	Marcon, G., Naveau, P. and Padoan, S. A. (2017).
	A semi-parametric stochastic generator for bivariate extreme events,
	\emph{Stat 6}(1), 184--201.
}

\examples{
################################################
# The following examples provide the left panels
# of Figure 1, 2 & 3 of Marcon et al. (2017).
################################################

## Figure 1 - symmetric logistic

\donttest{
# Strong dependence
a <- angular(model='log', n=50, dep=0.3, seed=4321, k=20, nsim=10000)
# Mild dependence
b <- angular(model='log', n=50, dep=0.6, seed=212, k=10, nsim=10000)
# Weak dependence
c <- angular(model='log', n=50, dep=0.9, seed=4334, k=6, nsim=10000)
}

## Figure 2 - Asymmetric logistic

\donttest{
# Strong dependence
d <- angular(model='alog', n=25, dep=0.3, asy=c(.3,.8), seed=43121465, k=20, nsim=10000)
# Mild dependence
e <- angular(model='alog', n=25, dep=0.6, asy=c(.3,.8), seed=1890, k=10, nsim=10000)
# Weak dependence
f <- angular(model='alog', n=25, dep=0.9, asy=c(.3,.8), seed=2043, k=5, nsim=10000)
}

}
