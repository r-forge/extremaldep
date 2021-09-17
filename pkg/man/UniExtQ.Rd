\name{UniExtQ}
\alias{UniExtQ}

\title{Univariate Extreme Quantile}
\description{Computes the extreme-quantiles of a univariate random variable corresponding to some exceedance probabilities.}
\usage{
UniExtQ(data, P=NULL, method="bayesian", U=NULL,
		cov=as.matrix(rep(1,length(data))), QatCov=NULL,
		param0=NULL, sig0=NULL, nsim=NULL, p=0.234,
		optim.meth="BFGS", control=NULL, ...)
}

\arguments{
	\item{data}{A vector of \eqn{n} observations.}
	\item{P}{The vector of probabilities associated to the quantiles.}
	\item{method}{The estimation method can be "bayesian" or "frequentist".}
	\item{U}{The threshold value under which the observations are censored.}
	\item{cov}{A \eqn{n \times c} matrix of covariates for the location parameter.}
	\item{QatCov}{q \eqn{n \times c} matrix with the value of the covariates at which the quantiles should be computed.}
	\item{param0}{The vector of initial value for the parameters. It should be of length \eqn{c + 2}}
	\item{sig0}{Only required when \code{method="bayesian"}. Initial value for the standard deviation of the multivarial normal proposal distribution.}
	\item{nsim}{Only required when \code{method="bayesian"}. Number of iterations in the Metropolis-Hastings algorithm.}
	\item{p}{Only required when \code{method="bayesian"}. Targeted acceptance ratio.}
	\item{optim.meth}{Only required when \code{method="frequentist"}. Optimisation algorithm as defined in \code{optim} function.}
	\item{control}{Only required when \code{method="frequentist"}. See details of \code{optim} function.}
	\item{...}{Only required when \code{method="frequentist"}. See details of \code{optim} function.}

}

\details{
	For some dataset given by \code{data}, the extreme-quantiles corresponding to some exceedance probability(ies) given in \code{P} are computed. The observations below the threshold \code{U} are considered censored.
	\itemize{
		\item If \code{cov} is specified then a linear model for the location parameter is fitted and \code{QatCov} must be specified. It sets the value of the covariates at which the extreme quantiles are evaluated;
		\item If \code{method=="frequentist"} then the GEV parameters are estimated via optimisation of the censored likelihood;
		\item If \code{method=="bayesian"} the the GEV parameters are estimated via a Metropolis-Hastings algorithm, where the proposal distirbution is a multivariate normal. After running for \code{nsim} iteration the algorithm is paused and some diagnostics plots are printed. The user then needs to enter the value of the burn-in period.
	}
	Refer to Beranger et al. (2019) for more details about the estimation procedures.
}

\value{
	If \code{method=="frequentist"}, a list with elements:
	\itemize{
		\item \code{Q.est}: a matrix of extreme quantiles where each of the \code{length(P)} rows represents an associated probability \code{P} and each of the \code{nrow(QatCov)} columns represents a level of the covariates;
		\item \code{kn}: the proportion of observation above the threshold \code{U}, corresponds to \eqn{\frac{k}{n}}{\frac{k}{n}};
		\item \code{est}: the vector of estimated parameters from the maximisation of the censored likelihood (on the log scale);
		\item \code{VarCov}: the variance-covariance matrix of the parameter estimates. Available if the \code{hessian=TRUE} is specified.
	}
	If \code{method=="bayesian"}, a list with elements:
	\itemize{
		\item \code{Q.est}: a list of extreme quantiles where each element is a vector and refers to a probability specified by \code{P};
		\item \code{post_sample}: a matrix containing the posterior sample of each parameter. The number of columns should be equal to \code{ncol(cov)+2};
		\item \code{straight.reject}: the number of straight rejections from the proposal distribution;
		\item \code{sig.vec}: the vector of updated scale parameter in the multivariate normal proposal.
	}
}

\references{
  Beranger, B., Padoan S. A. and Sisson, S. A. (2019).
  Estimation and uncertainty quantification for extreme quantile regions.
  \emph{arXiv e-prints} arXiv:1904:08251.
}


\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://mypage.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com/};
}

\seealso{
	\link{ExtQset}
}

\examples{
if (interactive()){

distribution <- "FRE" # MDA(FRE), Tail index "xi"

cat("\n Distribution:", distribution, "\n")
set.seed(999)
n <- 1500
loc <- 3
scale <- 1
shape <- 1/3
prob <- 0.9
cov <- as.matrix(rep(1,n)) # No covariates

data <- rfrechet(n = n, mu = loc, sigma = scale, lambda =shape)
pars <- c(loc, scale, shape)
pars.name <- paste("loc", loc*10, "_scale", scale*10, "_shape", shape*10, sep="")

### Required for Bayesian Estimation
nsim <- 5e+4
param0 <- c(1, 2, 1)
P <-c(1/750, 1/1500, 1/3000)
U <- quantile(data, probs=prob, type=3)
sig0 <- 1

### Estimation

# Bayesian approach:
mcmc.name <- paste("mcmc",distribution,"_n", n, "_prob", prob, "_", pars.name,sep="")
op <- par(mar=c(4.2,4.6,0.3,0.1))

assign(mcmc.name, UniExtQ(data=data, P=P, U=U, cov=cov, param0=param0, sig0=sig0,
	   nsim=nsim))
### RUN UNTIL HERE AND SPECIFY BURN-IN PERIOD

# Frequentist approach:
q <- UniExtQ(data=data, method="frequentist", P=P, param0=param0,
			 control = list(maxit = 5e+5, trace = 2))

### Illustration

ti <- 1/shape
mcmc <- get(mcmc.name)
Kern <- density(mcmc$post_sample[,ncol(cov)+2]) # KDE of the tail index

Hist <- hist(mcmc$post_sample[,ncol(cov)+2], prob=TRUE, col="lightgrey",
			 ylim=range(Kern$y), main="", xlab="Tail Index",
			 cex.lab=1.8, cex.axis=1.8, lwd=2)
ti_ic <- quantile(mcmc$post_sample[,ncol(cov)+2], probs=c(0.025, 0.975))
points(x=ti_ic, y=c(0,0), pch=4, lwd=4)
lines(Kern, lwd = 2, col = "dimgrey")
abline(v=ti, lwd=2)
abline(v=mean(mcmc$post_sample[,ncol(cov)+2]), lwd=2, lty=2)

#### Check the ability to estimate quantile regions
Q <- qfrechet(p = P, mu = loc, sigma = scale, lambda = shape, lower.tail=FALSE)

ci <- apply(log(mcmc$Q.est),2,function(x) quantile(x, probs=c(0.025, 0.975)))
for(i in 1:length(P)){
 	cat("\n Confidence interval extreme quantile with Probability ", P[i],
 		" : (", ci[1,i],",",ci[2,i],") \n")
}

Kern.est <- apply(log(mcmc$Q.est),2,density)

R <- range(log(c(Q,mcmc$Q.est,data)))
Xlim <- c(log(quantile(data,0.95)), R[2])
Ylim <- c(0, max(Kern.est[[1]]$y, Kern.est[[2]]$y, Kern.est[[3]]$y))

plot(log(data), rep(0,n), pch=16, main="", xlim=Xlim, ylim=Ylim,
	 xlab=expression(log(x)), ylab="Density", cex.lab=1.8, cex.axis=1.8, lwd=2)
polygon(Kern.est[[1]], col= rgb(211,211,211, 0.8*255, maxColorValue = 255),
		border=rgb(211,211,211, 255, maxColorValue = 255), lwd=4)
polygon(Kern.est[[2]], col= rgb(169,169,169, 0.8*255, maxColorValue = 255),
		border=rgb(169,169,169, maxColorValue = 255), lwd=4)
polygon(Kern.est[[3]], col= rgb(105,105,105, 0.8*255, maxColorValue = 255),
		border=rgb(105,105,105, maxColorValue = 255), lwd=4)
points(log(data), rep(0,n), pch=16, lwd=2)
abline(v=log(Q), lwd=2, lty=1)
for(j in 1:length(P)){abline(v=mean(log(mcmc$Q.est[,j])), lwd=2, lty=2)}
abline(v=log(q$Q.est), lwd=2, lty=3)

par(op)

}
}


\keyword{htest}
