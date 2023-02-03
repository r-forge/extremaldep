\name{ExtQ}
\alias{ExtQ}

\title{Univariate Extreme Quantile}
\description{Computes the extreme-quantiles of a univariate random variable corresponding to some exceedance probabilities.}
\usage{
ExtQ(P=NULL, method="Frequentist", pU=NULL, 
     cov=NULL, param=NULL,  param_post=NULL)
}


\arguments{
	\item{P}{A vector with values in \eqn{[0,1]}{[0,1]} indicating the probabilities of the quantiles to be computed.}
	\item{method}{A character string indicating the estimation method. Takes value \code{"bayesian"} or \code{"frequentist"}.}
	\item{pU}{A value in \eqn{[0,1]}{[0,1]} indicating the probability of exceeding a high threshold. In the estimation procedure, observations below the threshold are censored.}
	\item{cov}{A \eqn{q \times c}{q \times c} matrix indicating \eqn{q}{q} observations of \eqn{c-1}{c-1} covariates for the location parameter.}
	\item{param}{A \eqn{(c + 2)}{(c+2)} vector indicating the estimated parameters. Required when \code{method="Frequentist"}.}
	\item{param_post}{A \eqn{n \times (c + 2)}{n \times (c+2)} matrix indicating the posterior sample for the parameters, where \eqn{n}{n} is the number of MCMC replicates after removal of the burn-in period. Required when \code{method="Bayesian"}.}
}

\details{
  The first column of \code{cov} is a vector of 1s corresponding to the intercept. 

  When \code{pU} is \code{NULL} (default), then it is assumed that a block maxima approach was taken and quantiles are computed using the \code{\link{qGEV}} function. When \code{pU} is provided, the it is assumed that a threshold exceedances approach is taken and the quantiles are computed as
  \deqn{\mu + \sigma * \left(\left(\frac{pU}{P}\right)^\xi-1\right) \frac{1}{\xi}.}{\mu + \sigma * ((pU/P)^\xi-1) / \xi.}
}

\value{
	When \code{method=="frequentist"}, the function returns a vector of length \code{length(P)} if \code{ncol(cov)=1} (constant mean) or a \code{(length(P) x nrow(cov))} matrix if \code{ncol(cov)>1}. \cr
	
  When \code{method=="bayesian"}, the function returs a \code{(length(param_post) x length(P))} matrix if \code{ncol(cov)=1} or a list of \code{ncol(cov)} elements each taking a \code{(length(param_post) x length(P))} matrix if \code{ncol(cov)>1}.

}

\references{
	Beranger, B., Padoan, S. A. and Sisson, S. A. (2021).
	Estimation and uncertainty quantification for extreme quantile regions.
  \emph{Extremes}, \bold{24}, 349-375.
}


\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com}
}

\seealso{
		\code{\link{fGEV}}, \code{\link{qGEV}}
}

\examples{

##################################################
### Example - Pollution levels in Milan, Italy ###
##################################################

if(interactive()){

data(MilanPollution)

# Frequentist estimation
fit <- fGEV(Milan.winter$PM10)
fit$est

q1 <- ExtQ(P=1/c(600,1200,2400), method="Frequentist", param=fit$est)
q1

# Bayesian estimation with high threshold
cov <- cbind(rep(1,nrow(Milan.winter)), Milan.winter$MaxTemp, 
             Milan.winter$MaxTemp^2)
u <- quantile(Milan.winter$PM10, prob=0.9, type=3, na.rm=TRUE)

fit2 <- fGEV(data=Milan.winter$PM10, par.start=c(50,0,0,20,1), 
               method="Bayesian", u=u, cov=cov, sig0=0.1, nsim=5e+4) 

r <- range(Milan.winter$MaxTemp, na.rm=TRUE)
t <- seq(from=r[1], to=r[2], length=50)
pU <- mean(Milan.winter$PM10>u, na.rm=TRUE)
q2 <- ExtQ(P=1/c(600,1200,2400), method="Bayesian", pU=pU,
           cov=cbind(rep(1,50), t, t^2),
           param_post=fit2$param_post[-c(1:3e+4),])
             
R <- c(min(unlist(q2)), 800)
qseq <- seq(from=R[1],to=R[2], length=512)
Xl <- "Max Temperature"
Yl <- expression(PM[10])
  
for(i in 1:length(q2)){
  K_q2 <- apply(q2[[i]],2, function(x) density(x, from=R[1], to=R[2])$y)
  D <- cbind(expand.grid(t, qseq), as.vector(t(K_q2)) )
  colnames(D) <- c("x","y","z")
  fields::image.plot(x=t, y=qseq, z=matrix(D$z, 50, 512), xlim=r, 
                           ylim=R, xlab=Xl, ylab=Yl)
}

}
  
##########################################################
### Example - Simulated data from Frechet distirbution ###
##########################################################
  
if(interactive()){  
  
set.seed(999)  
data <- extraDistr::rfrechet(n=1500, mu=3, sigma=1, lambda=1/3)

u <- quantile(data, probs=0.9, type=3)
fit3 <- fGEV(data=data, par.start=c(1,2,1), method="Bayesian", 
             u=u, sig0=1, nsim=5e+4)
  
pU <- mean(data>u)
P <- 1/c(750,1500,3000)
q3 <- ExtQ(P=P, method="Bayesian", pU=pU,
           param_post=fit3$param_post[-c(1:3e+4),])  

### Illustration

# Tail index estimation

ti_true <- 3
ti_ps <- fit3$param_post[-c(1:3e+4),3]

K_ti <- density(ti_ps) # KDE of the tail index
H_ti <- hist(ti_ps, prob=TRUE, col="lightgrey",
			       ylim=range(K_ti$y), main="", xlab="Tail Index",
			       cex.lab=1.8, cex.axis=1.8, lwd=2)
ti_ic <- quantile(ti_ps, probs=c(0.025, 0.975))

points(x=ti_ic, y=c(0,0), pch=4, lwd=4)
lines(K_ti, lwd = 2, col = "dimgrey")
abline(v=ti_true, lwd=2)
abline(v=mean(ti_ps), lwd=2, lty=2)  
  
# Quantile estimation

q3_true <- extraDistr::qfrechet(p=P, mu=3, sigma=1, lambda=1/3, lower.tail=FALSE)

ci <- apply(log(q3), 2, function(x) quantile(x, probs=c(0.025, 0.975)))
K_q3 <- apply(log(q3), 2, density)

R <- range(log(c(q3_true, q3, data)))
Xlim <- c(log(quantile(data, 0.95)), R[2])
Ylim <- c(0, max(K_q3[[1]]$y, K_q3[[2]]$y, K_q3[[3]]$y))

plot(0, main="", xlim=Xlim, ylim=Ylim, xlab=expression(log(x)), 
     ylab="Density", cex.lab=1.8, cex.axis=1.8, lwd=2)
cval <- c(211, 169, 105)	 
for(j in 1:length(P)){
  col <- rgb(cval[j], cval[j], cval[j], 0.8*255, maxColorValue=255)
  col2 <- rgb(cval[j], cval[j], cval[j],  255, maxColorValue=255)
  polygon(K_q3[[j]], col=col, border=col2, lwd=4)
}
points(log(data), rep(0,n), pch=16)	 
# add posterior means
abline(v=apply(log(q3),2,mean), lwd=2, col=2:4)
# add credible intervals
abline(v=ci[1,], lwd=2, lty=3, col=2:4)
abline(v=ci[2,], lwd=2, lty=3, col=2:4)

}

}

\keyword{htest}
