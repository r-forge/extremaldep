\name{AngDensPlot}
\alias{AngDensPlot}

\title{3-D plot of parametric angular densities.}
\description{Plots (log)-angular densities on the three-dimensional simplex. Contour levels and data points (optional) are represented.}

\usage{
AngDensPlot(model='Pairwise', para=c(2,4,15,1), log=TRUE, data=NULL,
            contour=TRUE, labels=c(expression(w[1]),expression(w[3]),
            expression(w[2])), cex.dat=1, cex.lab=1, cex.cont=1)
}

\arguments{
	\item{model}{A string with the name of the parametric model for the angular density. \code{"Pairiwse"} is the default,
  see \code{\link{alik}}.}
	\item{para}{A numeric vector with the parameters of the parameteric model. Default is \code{para=c(2,4,15,1)},
		parameters for the Pairwise Beta model.}
	\item{log}{Logical; if \code{log=TRUE} (default) then the log-density is plotted.}
	\item{data}{If a (three-dimensional) dataset if provided then the data points are added to the density plot.}
	\item{contour}{Logical; if \code{contour=TRUE} (default) then the contour levels are plotted.}
	\item{labels}{Labels for the three corners of the simplex. Default is
		\code{labels=c(expression(w[1]),expression(w[3]),expression(w[2]))}. See \bold{details}.}
	\item{cex.dat}{Magnification of data points. Only if \code{data != NULL}.}
	\item{cex.lab}{Magnification of the labels. }
	\item{cex.cont}{Magnification of the contour labels.}
}

\details{
	Contour levels are given for the deciles. If \code{data != NULL} then the deciles are calculated using the density
	values of each data point. If \code{data = NULL} then the deciles are calcuated using all the points of the grid.

	\code{labels} are given in an anticlockwise order: bottom right, top middle and bottom left.
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.unibocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{http://www.borisberanger.com};
}

\references{
	Beranger, B. and Padoan, S. A. (2015).
	Extreme dependence models, chapater of the book \emph{Extreme Value Modeling and Risk Analysis: Methods and Applications},
  \bold{Chapman Hall/CRC}.

  	Beranger, B., Padoan, S. A. and Sisson, S. A. (2017).
	Models for extremal dependence derived from skew-symmetric families.
	\emph{Scandinavian Journal of Statistics}, \bold{44}(1), 21-45.

}

\examples{
################################################
# The following examples provide the plots of
# Figure 1.2 of the paper Beranger and Padoan (2015)
################################################

# The code has been frozen to speed up the package check.
# Please remove the hash symbol to test the code.

# Asymmetric Logistic

\donttest{
AngDensPlot('Asymmetric', c(rep(1,3),5.75, rep(0,6), 0.5,0.5,0.5), FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot('Asymmetric', c(rep(1,3),1.01, rep(0,6), 0.9,0.9,0.9), FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot('Asymmetric', c(rep(1,3),1.25, rep(0,6), 0.5,0.5,0.5), FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot('Asymmetric', c(rep(1,3),1.4, rep(0,6),0.7,0.15,0.15), FALSE, cex.lab=1.5, cex.cont=1.3)
}

# Tilted Dirichlet

\donttest{
AngDensPlot(model='Dirichlet', para=c(2,2,2), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Dirichlet', para=c(0.5,0.5,0.5), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Dirichlet', para=c(2,2.5,30), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Dirichlet', para=c(0.1,0.25,0.95), log=FALSE, cex.lab=1.5, cex.cont=1.3)
}

# Pairwise Beta

\donttest{
AngDensPlot(model='Pairwise', para=c(2,2,2,4), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Pairwise', para=c(1,1,1,0.5), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Pairwise', para=c(2,4,15,1), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Pairwise', para=c(10,10,10,1), log=FALSE, cex.lab=1.5, cex.cont=1.3)
}

# Husler-Reiss

\donttest{
AngDensPlot(model='Husler', para=c(0.3,0.3,0.3), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Husler', para=c(1.4,1.4,1.4), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Husler', para=c(1.7,0.7,1.1), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Husler', para=c(0.52,0.71,0.52), log=FALSE, cex.lab=1.5, cex.cont=1.3)
}

# Extremal-t

\donttest{
AngDensPlot(model='Extremalt', para=c(0.95,0.95,0.95,2), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Extremalt', para=c(-0.3,-0.3,-0.3,5), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Extremalt', para=c(0.52,0.71,0.52,3), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Extremalt', para=c(0.52,0.71,0.52,2), log=FALSE, cex.lab=1.5, cex.cont=1.3)
}

################################################
# The following examples provide
# the plots of Figure 1.3 of the paper
# Beranger and Padoan (2015)
################################################

## Load datasets
data(pollution)

Nsim <- 50e+4
Nbin <- 30e+4
MCpar <- 0.35
Hpar.pb <- 	list(mean.alpha=0, mean.beta=3,sd.alpha=3, sd.beta=3)
Hpar.hr <- list(mean.lambda=0, sd.lambda=3)
Hpar.di <- list(mean.alpha=0, sd.alpha=3)
Hpar.et <- list(mean.rho=0, mean.mu=3,sd.rho=3, sd.mu=3)

## PNS Data

\donttest{
est.pb.PNS <- posteriorMCMC(Nsim, Nbin, Hpar.pb, MCpar, PNS, seed=14342, model='Pairwise')
est.hr.PNS <- posteriorMCMC(Nsim, Nbin, Hpar.hr, MCpar, PNS, seed=14342, model='Husler')
est.di.PNS <- posteriorMCMC(Nsim, Nbin, Hpar.di, MCpar, PNS, seed=14342, model='Dirichlet')

lab1 <- c("PM10","NO","SO2")
AngDensPlot("Pairwise", est.pb.PNS$emp.mean, data=PNS, labels=lab1, cex.dat=0.8)
AngDensPlot("Husler", est.hr.PNS$emp.mean, data=PNS, labels=lab1, cex.dat=0.8)
AngDensPlot("Dirichlet", est.di.PNS$emp.mean, data=PNS, labels=lab1, cex.dat=0.8)
}

## NSN data

\donttest{
est.pb.NSN <- posteriorMCMC(Nsim, Nbin, Hpar.pb, MCpar, NSN, seed=14342, model='Pairwise')
est.hr.NSN <- posteriorMCMC(Nsim, Nbin, Hpar.hr, MCpar, NSN, seed=14342, model='Husler')
est.di.NSN <- posteriorMCMC(Nsim, Nbin, Hpar.di, MCpar, NSN, seed=14342, model='Dirichlet')

lab2 <- c("NO2","NO","SO2")
AngDensPlot("Pairwise", est.pb.NSN$emp.mean, data=NSN, labels=lab2, cex.dat=0.8)
AngDensPlot("Husler", est.hr.NSN$emp.mean, data=NSN, labels=lab2, cex.dat=0.8)
AngDensPlot("Dirichlet", est.di.NSN$emp.mean, data=NSN, labels=lab2, cex.dat=0.8)
}

## PNN data

\donttest{
est.pb.PNN <- posteriorMCMC(Nsim, Nbin, Hpar.pb, MCpar, PNN, seed=14342, model='Pairwise')
est.hr.PNN <- posteriorMCMC(Nsim, Nbin, Hpar.hr, MCpar, PNN, seed=14342, model='Husler')
est.di.PNN <- posteriorMCMC(Nsim, Nbin, Hpar.di, MCpar, PNN, seed=14342, model='Dirichlet')

lab3 <- c("PM10","NO","NO2")
AngDensPlot("Pairwise", est.pb.PNN$emp.mean, data=PNN, labels=lab3, cex.dat=0.8)
AngDensPlot("Husler", est.hr.PNN$emp.mean, data=PNN, labels=lab3, cex.dat=0.8)
AngDensPlot("Dirichlet", est.di.PNN$emp.mean, data=PNN, labels=lab3, cex.dat=0.8)
}

################################################
# The following examples provide the plots of
# the supplementary material for
# Beranger, Padoan and Sisson (2016)
################################################

\donttest{
AngDensPlot(model='Skewt', para=c(0.6,0.8,0.7,0,0,0,3), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Skewt', para=c(0.6,0.8,0.7,-3,-3,7,3), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Skewt', para=c(0.6,0.8,0.7,7,-10,3,3), log=FALSE, cex.lab=1.5, cex.cont=1.3)

AngDensPlot(model='Skewt', para=c(0.7,0.7,0.7,0,0,0,3), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Skewt', para=c(0.7,0.7,0.7,-3,-3,7,3), log=FALSE, cex.lab=1.5, cex.cont=1.3)
AngDensPlot(model='Skewt', para=c(0.7,0.7,0.7,7,-10,3,3), log=FALSE, cex.lab=1.5, cex.cont=1.3)
}

}

\keyword{hplot}
