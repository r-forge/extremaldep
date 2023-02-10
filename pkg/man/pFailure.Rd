\name{pFailure}
\alias{pFailure}

\title{Probability of falling into a failure region}
\description{This function computes the empirical estimate of the probability of falling into two types of failure regions.}

\usage{
pFailure(n, beta, u1, u2, mar1, mar2, type, plot, xlab, ylab, 
         nlevels=10)
}

\arguments{
	\item{n}{An integer indicating the number of points generated to compute the empirical probability.}
	\item{beta}{A vector representing the Bernstein polynomial coefficients.}
	\item{u1, u2}{Vectors of positive reals representing some high thresholds.}
	\item{mar1, mar2}{Vectors of marginal (GEV) parameters}
	\item{type}{A character string indicating if the failure region includes at least one exceedance (\code{"or"}), or both exceednaces (\code{"and"}). If \code{"both"} then probabilities to fall in both regions are computed.}
	\item{plot}{A logical value; if \code{TRUE} contour plots of the probalities are displayed.}
	\item{xlab, ylab}{A character string equivalent representing the graphical parameters as in \code{par}.}
	\item{nlevels}{The number of contour levels to be displayed.}
}

\details{
The two failure regions are:
\deqn{A_u = \left\{ (v_1, v_2): v_1>u_1 \textrm{ or } v_2 > u_2\right\} }{Au = {(v1,v2): v1>u1 or v2>u2}}
and
\deqn{B_u = \left\{ (v_1, v_2): v_1>u_1 \textrm{ and } v_2 > u_2\right\} }{Bu = {(v1,v2): v1>u1 and v2>u2}}
where \eqn{(v_1,v_2) \in R_+^2}{(v1,v2) \in R_+^2} and \eqn{u_1, u_2>0}{u1, u2>0}.

Exceedances samples \eqn{y_{1,1}, \ldots, y_{n_1}}{y11, ..., yn1} and \eqn{y_{1,2}, \ldots, y_{n_2}}{y12, ..., yn2} are generating according to Algorithm 3 of Marcon et al. (2017) and the the estimates of the probability of falling in \eqn{A_u}{Au} and \eqn{B_u}{Bu} are given by
\deqn{\hat{P}_{A_u} = \frac{1}{n}\sum{i=1}^n I\left\{ y_{i,1}> u_1^* \textrm{ or } y_{i,2}> u_2^* \right\} }{P^Au = 1/n * \sum_{i=1}^n I{yi1 >u1*  or yi2 > u2*}}
and
\deqn{\hat{P}_{B_u} = \frac{1}{n}\sum{i=1}^n I\left\{ y_{i,1}> u_1^* \textrm{ and } y_{i,2}> u_2^* \right\} }{P^Bu = 1/n * \sum_{i=1}^n I{yi1 >u1* and yi2 > u2*}.}
}

\value{
	A list containing \code{AND} and/or \code{OR}, depending on the \code{type} argument. Each element is a \code{length(u1) x length(u2)} matrix.   
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}

\references{

	Marcon, G., Naveau, P. and Padoan, S.A. (2017)
  A semi-parametric stochastic generator for bivariate extreme events
  \emph{Stat}, \bold{6}, 184-201.
}

\seealso{\code{\link{dExtDep}}, \code{\link{rExtDep}}, \code{\link{fExtDep}}, \code{\link{fExtDep.np}} }

\examples{

# Example wind speed and wind gust

data(WindSpeedGust)
years <- format(ParcayMeslay$time, format="\%Y")
attach(ParcayMeslay[which(years \%in\% c(2004:2013)),])
  
WS_th <- quantile(WS,.9)
DP_th <- quantile(DP,.9)
  
pars.WS <- evd::fpot(WS, WS_th, model="pp")$estimate
pars.DP <- evd::fpot(DP, DP_th, model="pp")$estimate
  
data_uf <- trans2UFrechet(cbind(WS,DP), type="Empirical")
  
rdata <- rowSums(data_uf)
r0 <- quantile(rdata, probs=.90)
extdata <- data_uf[rdata>=r0,]
  
SP_mle <- fExtDep.np(method="Frequentist", data=extdata, k0=10, 
                     type="maxima")

if(interactive()){
  pF <- pFailure(n=50000, beta=SP_mle$Ahat$beta,
                 u1=seq(from=19, to=28, length=200), mar1=pars.WS,
                 u2=seq(from=40, to=60, length=200), mar2=pars.DP, 
                 type="both", plot=TRUE, 
                 xlab="Daily-maximum Wind Speed (m/s)",
                 ylab="Differential of Pressure (mbar)", nlevels=15)
}

}

\keyword{models}
