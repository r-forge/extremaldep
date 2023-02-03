\name{fExtDepSpat}
\alias{fExtDepSpat}

\title{Fitting of a max-stable process}
\description{This function uses the Stephenson-Tawn likelihood to estimate parameters of max-stable models.}

\usage{
fExtDepSpat(model, z, sites, hit, jw, thresh, DoF, range, smooth, 
            alpha, par0, acov1, acov2, parallel, ncores, args1, args2, 
            seed=123, method = "BFGS", sandwich=TRUE,
            control = list(trace=1, maxit=50, REPORT=1, reltol=0.0001))
}

\arguments{
	\item{model}{A character string indicating the max-stable model, currently extremal-t (\code{"ET"}) and extremal skew-t (\code{"EST"}) only available. Note that the Schlather model can be obtained as a special case by specifying the extremal-t model with \code{DoF=1}}
	\item{z}{A \eqn{(n \times d)}{(n x d)} matrix containing \eqn{n}{n} observations at \eqn{d}{d} locations.}
	\item{sites}{A \eqn{(d \times 2)}{(d x 2)} matrix corresponding to the coordinates of locations where the processes is simulated. Each row corresponds to a location.}
	\item{hit}{A \eqn{(n \times d)}{(n x d)} matrix containing the hitting scenarios for each observations.}
	\item{jw}{An integer between \eqn{2}{2} and \eqn{d}{d}, indicating the tuples considered in the composite likelihood. If \code{jw=d} the full likelihood is considered.}
	\item{thresh}{A positive real indicating the threshold value for pairwise distances. See \code{details}.}
	\item{DoF}{A positive real indicating a fixed value of the degree of freedom of the extremal-t and extremal skew-t models.}
	\item{range}{A positive real indicating a fixed value of the range parameter for the power exponential correlation function (only correlation function currently available).}
	\item{smooth}{A positive real in \eqn{(0,2]}{(0,2])} indicating a fixed value of the smoothness parameter for the power exponential correlation function (only correlation function currently available).}
	\item{alpha}{A vector of length \eqn{3}{3)} indicating fixed values of the skewness parameters \eqn{\alpha_0}{\alpha_0}, \eqn{\alpha_1}{\alpha_1} and \eqn{\alpha_2}{\alpha_2} for the extremal skew-t model. If some components are \code{NA}, then the corresponding parameter will be estimated.}
	\item{par0}{A vector of initial value of the parameter vector, in order the degree of freedom \eqn{\nu}{\nu}, the range \eqn{r}{r}, the smoothness \eqn{\eta}{\eta} and the skewness parameters \eqn{\alpha_0}{\alpha_0}, \eqn{\alpha_1}{\alpha_1}. Its length depends on the model and the number of fixed parameters.}
	\item{acov1, acov2}{Vectors of length \eqn{d}{d} representing covariates to model the skewness parameter of the extremal skew-t model.}
	\item{parallel}{A logical value; if \code{TRUE} parallel computing is done.}
	\item{ncores}{An integer indicating the number of cores considered in the parallel socket cluster of type \code{'PSOCK'}, based on the \code{makeCluster} routine from the \code{parallel} package. Required if \code{parallel=TRUE}.}
	\item{args1, args2}{Lists specifying details about the Monte Carlo simulation schereme to compute multivariate CDFs. See \code{details}.}
	\item{seed}{An integer for reproduciblity in the CDF computations.}
	\item{method}{A character string indicating the optimisation method to be used. See \code{\link{optim}} for more details.}
	\item{sandwich}{A logical value; if \code{TRUE} the standard errors of the estimates are computed from the Sandwich (Godambe) information matrix.}
	\item{control}{A list of control parameter for the optimisation. See \code{\link{optim}} for more details.}
}

\details{
  This routine follows the methodology developped by Beranger et al. (2021). It uses on the Stephenson-Tawn which relies on the knowledge of time occurrences of each block maxima. Rather than considering all partitions of the set \eqn{\{1, \ldots,d\}}{{1,\ldots,d}}, the likelihood is computed using the observed partition. Let \eqn{\Pi = (\pi_1, \ldots, \pi_K)}{\Pi = (\pi_1, \ldots, \pi_K)} denote the observed partition, then the Stephenson-Tawn likelihood is given by
\deqn{L(\theta;z) = \exp \left\{ - V(z;\theta) \right\} \times \prod^{K}_{k=1} - V_{\pi_k}(z;\theta),}{L(\theta;z) = \exp {- V(z;\theta)} x \prod^{K}_{k=1} - V_{\pi_k}(z;\theta),} where \eqn{V_{\pi}}{V_\pi} represents the partial derivative(s) of \eqn{V(z;\theta)}{V(z;\theta)} with respect to \eqn{\pi}{\pi}.

  When \code{jw=d} the full Stephenson-Tawn likelihood is considered whereas for values lower than \eqn{d}{d} a composite likelihood approach is taken. 
  
  The argument \code{thresh} is required when the composite likelihood is used. A tuple of size \code{jw}, is assigned a weight of one if the maximum pairwise distance between corresponding locations is less that \code{thresh} and a weight of zero otherwise.
  
  Arguments \code{args1} and \code{args2} relate to specifications of the Monte Carlo simulation scheme to compute multivariate CDF evaluations. These should take the form of lists including the minimum and maximum number of simulations used (\code{Nmin} and \code{Nmax}), the absolute error (\code{eps}) and whether the error should be controlled on the log-scale (\code{logeps}).
}

\value{
  A list comprising of the vector of estimated parameters (\code{est}), the composite likelihood order (\code{jw}), the maximised log-likelihood value (\code{LL}). In addition, if \code{sandwich=TRUE} the the standard errors from the sandwich information matrix are reported via \code{stderr.sand} as well as the TIC for model selection (\code{TIC}). Finally, if the composite likelihood is considered, a matrix with all tuples considered with a weight of 1 are reported in \code{cmat}. 
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
set.seed(14342)

# Simulation of 20 locations
Ns <- 20
sites <- matrix(runif(Ns*2)*10-5,nrow=Ns,ncol=2)
for(i in 1:2) sites[,i] <- sites[,i] - mean(sites[,i])

# Simulation of 50 replicates from the Extremal-t model
Ny <- 50
z <- rExtDepSpat(Ny, sites, model="ET", cov.mod="powexp", DoF=1, 
                 range=3, nugget=0, smooth=1.5, 
                 control=list(method="exact"))
                 
# Fit the extremal-t using the full Stephenson-Tawn likelihood
args1 <- list(Nmax=50L, Nmin=5L, eps=0.001, logeps=FALSE)
args2 <- list(Nmax=500L, Nmin=50L, eps=0.001, logeps=TRUE)
if(interactive()){
fit1 <- fExtDepSpat(model="ET", z=z$vals, sites=sites, hit=z$hits,
                    par0=c(3,1,1), parallel=TRUE, ncores=6,
                    args1=args1, args2=args2, control = list(trace=0))
fit1$est
}
}

\keyword{models}
