\name{AngDensPlot}
\alias{AngDensPlot}
\alias{print.AngDensPlot}

\title{3-D plot of parametric angular densities.}
\description{Plot (log)-angular densities on the three-dimensional simplex. Contour levels and data points (optional) are represented.}

\usage{
AngDensPlot(model='Pairwise', para=c(2,4,15,1), log=TRUE, data=NULL, 
            contour=TRUE, labels=c(expression(w[1]),expression(w[3]),
            expression(w[2])), cex.dat=1, cex.lab=1, cex.cont=1)
}

\arguments{
	\item{model}{A string with the name of the parametric model for the angular density. \code{"Pairiwse"} is the default, 
  see \code{\link{alik}}.}
	\item{para}{A numeric vector with the parameters of the parameteric model. Default is \code{para=c(2,4,15,1)},
		parameters for the Pairiwse Beta model.}	
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
	
	\code{labels} are given in the anticlockwise order: bottom right, top middle and bottom left. 
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.bocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
}

\references{  
	Beranger, B. and Padoan, S. A. (2015).
	Extreme dependence models, chapater of the book \emph{Extreme Value Modeling and Risk Analysis: Methods and Applications},
  \bold{Chapman Hall/CRC}.
}

\examples{
################################################
# The following esamples provide the plots of
# Figure 1.2 of the paper Beranger and Padoan (2015)
################################################

# The code has been frozen to speed up the package check. 
# Please remove the hash symbol to test the code.

# Asymmetric Logistic

# AngDensPlot(model='Asymmetric', para=c(5.75,0.5,0.5,0.5), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Asymmetric', para=c(1.01,0.9,0.9,0.9), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Asymmetric', para=c(1.25,0.5,0.5,0.5), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Asymmetric', para=c(1.4,0.7,0.15,0.15), log=FALSE, cex.lab=1.5, cex.cont=1.3)

# Tilted Dirichlet

# AngDensPlot(model='Dirichlet', para=c(2,2,2), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Dirichlet', para=c(0.5,0.5,0.5), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Dirichlet', para=c(2,2.5,30), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Dirichlet', para=c(0.1,0.25,0.95), log=FALSE, cex.lab=1.5, cex.cont=1.3)

# Pairwise Beta

# AngDensPlot(model='Pairwise', para=c(2,2,2,4), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Pairwise', para=c(1,1,1,0.5), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Pairwise', para=c(2,4,15,1), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Pairwise', para=c(10,10,10,1), log=FALSE, cex.lab=1.5, cex.cont=1.3)

# Husler-Reiss

# AngDensPlot(model='Husler', para=c(0.3,0.3,0.3), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Husler', para=c(1.4,1.4,1.4), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Husler', para=c(1.7,0.7,1.1), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Husler', para=c(0.52,0.71,0.52), log=FALSE, cex.lab=1.5, cex.cont=1.3)

# Extremal-t

# AngDensPlot(model='Extremalt', para=c(0.95,0.95,0.95,2), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Extremalt', para=c(-0.3,-0.3,-0.3,5), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Extremalt', para=c(0.52,0.71,0.52,3), log=FALSE, cex.lab=1.5, cex.cont=1.3)
# AngDensPlot(model='Extremalt', para=c(0.52,0.71,0.52,2), log=FALSE, cex.lab=1.5, cex.cont=1.3)

################################################
# The following esamples provide the plots of
# the plots of Figure 1.3 of the paper
# Beranger and Padoan (2015)
################################################

## Load datasets
# data(pollution) 

## PNS Data

# est.pb.PNS <- posteriorMCMC(Nsim=50e+4,Nbin=30e+4,dat= PNS, model='Pairwise',Hpar=list(mean.alpha=0, mean.beta=3,sd.alpha=3, sd.beta=3), MCpar=0.35,seed=14342)
# est.hr.PNS <- posteriorMCMC(Nsim=50e+4,Nbin=30e+4,dat= PNS, model='Husler',Hpar=list(mean.lambda=0, sd.lambda=3), MCpar=0.35,seed=14342)
# est.di.PNS <- posteriorMCMC(Nsim=50e+4,Nbin=30e+4,dat= PNS, model='Dirichlet',Hpar=list(mean.alpha=0,sd.alpha=3), MCpar=0.35,seed=14342)

# AngDensPlot(model="Pairwise", para=est.pb.PNS$emp.mean, log=TRUE, data=PNS, contour=TRUE,labels=c("PM10","NO","SO2"),cex.dat=0.8)
# AngDensPlot(model="Husler", para=est.hr.PNS$emp.mean, log=TRUE, data=PNS, contour=TRUE,labels=c("PM10","NO","SO2"),cex.dat=0.8)
# AngDensPlot(model="Dirichlet", para=est.di.PNS$emp.mean, log=TRUE, data=PNS, contour=TRUE,labels=c("PM10","NO","SO2"),cex.dat=0.8)

## NSN data

# est.pb.NSN <- posteriorMCMC(Nsim=50e+4,Nbin=30e+4,dat= NSN,Hpar=list(mean.alpha=0, mean.beta=3,sd.alpha=3, sd.beta=3), MCpar=0.35, model='Pairwise',seed=14342)
# est.hr.NSN <- posteriorMCMC(Nsim=50e+4,Nbin=30e+4,dat= NSN, model='Husler',Hpar=list(mean.lambda=0, sd.lambda=3), MCpar=0.35,seed=14342)
# est.di.NSN <- posteriorMCMC(Nsim=50e+4,Nbin=30e+4,dat= NSN, model='Dirichlet',Hpar=list(mean.alpha=0, sd.alpha=3), MCpar=0.35,seed=14342)

# AngDensPlot(model="Pairwise", para=est.pb.NSN, log=TRUE$emp.mean, data=NSN, contour=TRUE,labels=c("NO2","NO","SO2"),cex.dat=0.8)
# AngDensPlot(model="Husler", para=est.hr.NSN, log=TRUE$emp.mean, data=NSN, contour=TRUE,labels=c("NO2","NO","SO2"),cex.dat=0.8)
# AngDensPlot(model="Dirichlet", para=est.di.NSN, log=TRUE$emp.mean, data=NSN, contour=TRUE,labels=c("NO2","NO","SO2"),cex.dat=0.8)

## PNN data

# est.pb.PNN <- posteriorMCMC(Nsim=50e+4,Nbin=30e+4,dat= PNN,Hpar=list(mean.alpha=0, mean.beta=3,sd.alpha=3, sd.beta=3), MCpar=0.35, model='Pairwise',seed=14342)
# est.hr.PNN <- posteriorMCMC(Nsim=50e+4,Nbin=30e+4,dat= PNN, model='Husler',Hpar=list(mean.lambda=0, sd.lambda=3), MCpar=0.35,seed=14342)
# est.di.PNN <- posteriorMCMC(Nsim=50e+4,Nbin=30e+4,dat= PNN, model='Dirichlet',Hpar=list(mean.alpha=0, sd.alpha=3), MCpar=0.35,seed=14342)

# AngDensPlot(model="Pairwise", para=est.pb.PNN$emp.mean, log=TRUE, data=PNN, contour=TRUE,labels=c("PM10","NO","NO2"),cex.dat=0.8)
# AngDensPlot(model="Husler", para=est.hr.PNN$emp.mean, log=TRUE, data=PNN, contour=TRUE,labels=c("PM10","NO","NO2"),cex.dat=0.8)
# AngDensPlot(model="Dirichlet", para=est.di.PNN$emp.mean, log=TRUE, data=PNN, contour=TRUE,labels=c("PM10","NO","NO2"),cex.dat=0.8)
}