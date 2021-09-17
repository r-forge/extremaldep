\name{ExtQset}
\alias{ExtQset}

\title{Bivariate Extreme Quantile Sets}
\description{Computes extreme-quantiles regions of a bivariate random variable correspoding to some exceedance probabilities.}
\usage{
ExtQset(data, P=NULL, method="bayesian", U=NULL,
        cov1=as.matrix(rep(1,nrow(data))), cov2=as.matrix(rep(1,nrow(data))),
        QatCov1=NULL, QatCov2=NULL, mar=TRUE, par10=c(1,2,1), par20=c(1,2,1),
        sig10=1, sig20=1, param0=NULL, k0=NULL, pm0=NULL, prior.k="nbinom",
        prior.pm="unif", hyperparam = list(mu.nbinom = 3.2, var.nbinom = 4.48),
        nsim=NULL, lo=NULL, up=NULL, d=5)
}

\arguments{
	\item{data}{A matrix of \eqn{n \times 2} observations.}
	\item{P}{The vector of probabilities associated to the quantiles.}
	\item{method}{The estimation method can be "bayesian", "EdHK" or "frequentist".}
	\item{U}{The bivariate threshold value under which the observations are marginally censored.}
	\item{cov1}{A \eqn{n \times c1} matrix of covariates for the location parameter of the first margin.}
	\item{QatCov1}{q \eqn{n \times c1} matrix with the value of the first margin covariates at which the quantiles should be computed.}
	\item{cov2}{A \eqn{n \times c2} matrix of covariates for the location parameter of the second margin.}
	\item{QatCov2}{q \eqn{n \times c2} matrix with the value of the second margin covariates at which the quantiles should be computed.}
	\item{mar}{Only required when \code{method="bayesian"}. If \code{mar=TRUE} then a first estimation of the margins is done.}
	\item{par10, par20}{Only required when \code{method="bayesian"}. The vector of initial value for the parameters.}
	\item{sig10, sig20}{Only required when \code{method="bayesian"}. Initial value for the standard deviations of the multivariate normal proposal distribution for both margins.}
	\item{param0}{Only required when \code{method="bayesian"}. The vector of initial value for the Bernstein polynomial coefficients. It should be a list with elements \code{$eta} and \code{$beta}.}
	\item{k0}{Only required when \code{method="bayesian"}. The initial value of the polynomial order.}
	\item{pm0}{Only required when \code{method="bayesian"}. The list of initial values of the probability masses at the boundaries of the simplex. It should be a list with two elements \code{$p0} and \code{$p1}, see \link{bbeed}.}
	\item{prior.k}{Only required when \code{method="bayesian"}. The prior on the polynomial order, see \link{bbeed}.}
	\item{prior.pm}{Only required when \code{method="bayesian"}. The prior on the probability masses at the endpoints of the simplex, see \link{bbeed}.}
	\item{hyperparam}{Only required when \code{method="bayesian"}. A list of the hyper-parameters, see \link{bbeed}.}
	\item{nsim}{Only required when \code{method="bayesian"}. Number of iterations in the Metropolis-Hastings algorithm.}
	\item{lo}{Only required when \code{method="EdHK", "frequentist"}. Lower value of k in Hill estimator for shape parameter.}
	\item{up}{Only required when \code{method="EdHK", "frequentist"}. Upper value of k in Hill estimator for shape parameter.}
	\item{d}{postive integer, indicating the order of Bernstein polynomials}

}

\details{
	For some dataset given by \code{data}, the extreme-quantiles corresponding to some exceedance probability(ies) given in \code{P} are computed. The observations below the threshold \code{U} are considered censored.
	\itemize{
		\item If \code{method="bayesian"}, the methodologies given by \link{bbeed} and \link{UniExtQ} are combined. The algorithm is a Trans-dimensional MCMC scheme as described in Algorithm 1 of Beranger et al. (2019). If \code{mar=TRUE} then the function \link{UniExtQ} is preliminarily applied to the margins to select some starting values updating \code{par10, par20, sig10} and \code{sig20}. After running for \code{nsim} iteration the algorithm is paused and some diagnostics plots (on the margins and polynomial degree) are printed. The user then needs to enter the value of the burn-in period. See \link{bbeed} for details about the prior and hyperparameters for the dependence structure.
		\item If \code{method="EdHK"}, then the methodology developped by Einmahl et al. (2013) is applied. This is a non-parametric approach where the marginal parameters are estimated first.The Hill estimator is used to estimate the marginal tail indexes.
		\item If \code{method="frequentist"}, then similar to the \code{"EdHK"} estimator, the marginal parameters are estimated first (in the same way) and then the \link{beed} function is used to estimate the dependence structure.
	}
}

\value{
	If \code{method=="bayesian"}, a list with elements:
	\itemize{
		\item \code{Qset_P1, \ldots}: \code{length(P)} lists of \eqn{100}{100} \eqn{2}{2} by \eqn{3}{3} matrices corresponding to the extreme quantile regions associated with probability \code{P}, using \eqn{100}{100} equidistant points in the unit simplex, for \eqn{3}{3} different levels: \eqn{0.05}{0.05}-quantile, mean and \eqn{0.95}{0.95}-quantile;
		\item \code{Qset_P1_post, \ldots}: \code{length(P)} lists of \eqn{100}{100} \eqn{2}{2} by \eqn{npost}{npost} matrices corresponding to the posterior samples of size \code{npost} of the extreme quantile regions associated with probability \code{P};
		\item \code{ghat}: a \eqn{3}{3} by \eqn{100}{100} matrix giving the \eqn{0.05}{0.05}-quantile, mean and \eqn{0.95}{0.95}-quantile of the inverse angular density \eqn{q^*}{q^*} obtained from the posterior sample;
		\item \code{Shat}: a \eqn{3}{3} by \eqn{100}{100} matrix giving the \eqn{0.05}{0.05}-quantile, mean and \eqn{0.95}{0.95}-quantile of the basic set \eqn{\mathcal{B}}{\mathcal{B}} obtained from the posterior sample;
		\item \code{nuShat}: the \eqn{0.05}{0.05}-quantile, mean and \eqn{0.95}{0.95}-quantile of the estimate of the basic set size \eqn{\nu(\mathcal{B})}{\nu(\mathcal{B})};
		\item \code{burn}: the burn-in period, specified by the user after observing the diagnostic plots;
	}
	If \code{method=="EDhK"}, a list with elements:
	\itemize{
		\item \code{xn_hat1, \ldots}: \code{length(P)} vectors of lentgh \eqn{101}{101}giving the x-axis values of the estimated quantiles associated with probability \code{P}.
		\item \code{yn_hat1, \ldots}: \code{length(P)} vectors of lentgh \eqn{101}{101}giving the y-axis values of the estimated quantiles associated with probability \code{P}.
	}
	If \code{method=="frequentist"}, a list with elements:
	\itemize{
		\item \code{hhat}: a vector of length \eqn{100}{100} giving the estimated angular density at \eqn{100}{100} equally spaced points on the unit simplex;
		\item \code{ghat}: a vector of length \eqn{100}{100} giving the corresponding estimated \eqn{1/q^*}{1/q^*} function;
		\item \code{Shat}: a \eqn{100}{100} by \eqn{2}{2} matrix of the corresponding estimated basic set \eqn{\mathcal{B}}{\mathcal{B}};
		\item \code{nuShat}: a real giving the estimate of the basic set size \eqn{\nu(\mathcal{B})}{\nu(\mathcal{B})}
		\item \code{Qhat}: a \eqn{100 \times 2 \times \code{length(P)}}{100 \times 2 \times \code{length(P)}} list corresponding to \code{length(P)} \eqn{100 \times 2}{100 \times 2} matrices representing the estimated extreme quantile regions associated with probability \code{P};
		\item \code{gamhat}: a bivariate vector of the estimated marginal tail indices;
		\item \code{uhat}: a bivariate vector of the estimated location parameters.
	}
}

\references{
  Beranger, B., Padoan S. A. and Sisson, S. A. (2019).
  Estimation and uncertainty quantification for extreme quantile regions.
  \emph{arXiv e-prints} arXiv:1904:08251.

  Einmahl, J. H. J., de Haan, L. and Krajina, A. (2013).
  Estimating extreme bivariate quantile regions.
  \emph{Extremes}, \bold{16}, 121-145.
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://mypage.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com/};
	Andrea Krajina, \email{akrajina@gmail.com}
}

\seealso{
	\link{UniExtQ}, \link{bbeed}, \link{beed}
}

\examples{
if (interactive()){
library(mvtnorm)

distribution <- "Cauchy"

par10 <- par20 <- c(1,2,1) # Initial marginal parameter values
sig10 <- sig20 <- 1 # Initial scale values in MVN proposal
prior.k <- "nbinom" # Prior distribution for polynomial degree
k0 <- 5 # Degree of the polynomial
prior.pm <- "unif" # Prior distribution for the point masses
pm0 <- list(p0=0, p1=0)
# Vector of hyperparameters for prior distribution:
hyperparam <- list(mu.nbinom = 3.2, var.nbinom = 4.48, a.unif = 0, b.unif = 0.1)

###
### Data simulation
###

n <- 1500 # Sample size
P <- c(1/750, 1/1500, 1/3000) # Vector of probabilities for extreme quantiles
prob <- c(0.9, 0.9) # To use to evaluate thresholds

# Dependence structure;
rho <- 0
sigma <- matrix(c(1,rho,rho,1), ncol=2)
df <- 1

# Compute quantiles for the Cauchy:
ell1 <- ellipse(prob=1-P[1], pos=TRUE)
ell2 <- ellipse(prob=1-P[2], pos=TRUE)
ell3 <- ellipse(prob=1-P[3], pos=TRUE)

realx1 <- ell1[,1]; realy1 <- ell1[,2]
realx2 <- ell2[,1]; realy2 <- ell2[,2]
realx3 <- ell3[,1]; realy3 <- ell3[,2]

# Data simulation (Cauchy)
set.seed(999)
data <- rmvt(5*n, sigma=sigma, df=df)
data <- data[data[,1]>0 & data[,2]>0, ]
data <- data[1:n, ]

# Threshold
U <- c(quantile(data[,1], probs = prob[1], type=3), quantile(data[,2], probs = prob[2], type=3))

###
### Estimation
###

Q <- ExtQset(data=data, P=P, U=U, par10=par10, par20=par20, sig10=sig10, sig20=sig20, pm0=pm0,
             k0=k0, prior.k=prior.k, prior.pm=prior.pm, hyperparam=hyperparam, nsim=50000)

Q.EDhK <- ExtQset(data=data, P=P, method="EDhK", lo=50, up=300)

w <- seq(0.00001, .99999, length=100) # define grid
gfun <- ((w^2+(1-w)^2)^(-3/2))^(1/3) # Compute the true g function
xT <- gfun*w # x-axis of Basic set
yT <- gfun*(1-w) # y-axis of Basic set

###
### Graphical representation
###

op <- par(mfrow=c(2,3), mar=c(3, 3, 0.5, 0.2), mgp=c(2,.8,0))

# Plot 1: Density of Exponent measure

ylim.pl1 <- c(0,1.7)
plot(w, gfun, type="l", xlab="w", ylab=expression(1/q[symbol("\052")](w)), ylim=ylim.pl1)
polygon(c(w, rev(w)), c(Q$ghat[3,], rev(Q$ghat[1,])), col="gray")
lines(w, Q$ghat[2,],col="gray0", lwd=2, lty=3)
lines(w, gfun, lwd=2)

# Plot 2: Basic-set S

xlim.pl2 <-c(0,1.5); ylim.pl2 <- c(0,1.5)
plot(xT,yT, pch=19, col=1, type="l", xlim=xlim.pl2, ylim=ylim.pl2,
	 xlab=expression(x[1]), ylab=expression(x[2]))
polygon(c(Q$Shat[,1,3], rev(Q$Shat[,1,1])), c(Q$Shat[,2,3], rev(Q$Shat[,2,1])), col="gray")
points(Q$Shat[,,2], type="l", col="gray0", lwd=2, lty=3)
lines(xT,yT,lwd=2)

# Plot 3: Data + quantile regions

xlim.pl3 <- c(0, 3500); ylim.pl3 <- c(0, 3500)
plot(data, xlim=xlim.pl3, ylim=ylim.pl3, pch=19, xlab=expression(x[1]), ylab=expression(x[2]))
points(realx1,realy1, type="l", lwd=2, lty=1)
points(realx2,realy2, type="l", lwd=2, lty=1)
points(realx3,realy3, type="l", lwd=2, lty=1)
lines(Q$Qset_P1_CovNum_1[,,2], lty=3, col="gray0", lwd=2)
lines(Q$Qset_P2_CovNum_1[,,2], lty=3, col="gray0", lwd=2)
lines(Q$Qset_P3_CovNum_1[,,2], lty=3, col="gray0", lwd=2)

# Plot 4,5,6: Quantile region with probability 1/750, 1/1500, 1/3000

xlim.pl46 <- c(0,7400); ylim.pl46 <- c(0,7400)
for(j in 1:3){
  tmp.name <- paste("Qset_P",j,"_CovNum_1",sep="")
  tmp.quant <- Q[[tmp.name]]

  plot(data, xlim=xlim.pl46, ylim=ylim.pl46, type="n", pch=19,
  	xlab=expression(x[1]), ylab=expression(x[2]))
  polygon(c(tmp.quant[,1,3], rev(tmp.quant[,1,1])),
  	c(tmp.quant[,2,3], rev(tmp.quant[,2,1])), col="gray")
  points(get(paste("realx",j,sep="")), get(paste("realy",j,sep="")), type="l", lty=1, lwd=2)
  lines(tmp.quant[,,2], lty=3, col="gray0", lwd=2)
  lines(Q.EDhK[[paste("xn_hat",j,sep="")]], Q.EDhK[[paste("yn_hat",j,sep="")]], lty=2, lwd=2)
}
par(op)
}
}

\keyword{htest}
