\name{alik}
\alias{alik}

\title{Approximate likelihood estimation of extremal dependence models.}

\description{Estimates the parameters of extremal dependence models. It also provides
	standard errors and TIC.}

\usage{
alik(data, model, parastart, c=NULL, trace=0, sig=3)
}

\arguments{
	\item{data}{A (\eqn{n \times d}{n x d}) matrix of angular components, where the rows represent \eqn{n}{n} independent points
		in the \eqn{d}{d}-dimensional unit simplex.}
	\item{model}{A string with the name of the parametric model to be estimated. See \bold{Details}.}
	\item{parastart}{A vector containing the starting values of the model's parameters for the maximisation of
		the log-approximate likelihood. See \bold{Details}.}
	\item{c}{A real value in \eqn{[0,1]}{[0,1]}, providing the decision rule to allocate a data point to a subset of the simplex.
		Only required for the Extremal-t, Extremal Skew-t and Asymmetric Logistic models.}
	\item{trace}{Non-negative integer. See \code{\link{optim}} for details. \code{trace=0}
		is the default.}
	\item{sig}{Non-negative integer. Provides the number of decimal places for the returned object.
		\code{sig=3} is the default.}
}

\details{
The available parametric extremal dependence models are:
\itemize{
\item The \bold{Pairwise Beta}, called with \code{model="Pairwise"}. The number of parameters is
\code{choose(d,2)+1};
\item The \bold{Husler-Reiss}, called with \code{model="Husler"}. The number of parameters is
\code{choose(d,2)};
\item The \bold{Tilted Dirichlet}, called with \code{model="Dirichlet"}. The number of parameters is
\eqn{d}{d};
\item The \bold{Extremal-t}, called with \code{model="Extremalt"}. The number of parameters is \code{choose(d,2)+1};
\item The \bold{Extremal Skew-t}, called with \code{model="Skewt"}. The number of parameters is \code{choose(d,2)+d+1};
\item The \bold{Asymmetric Logistic}, that can be called with \code{model="Asymmetric"}. The number of dependence parameters is \eqn{2^{d-1}(d+2)-(2d+1)}{2^{d-1}(d+2)-(2d+1)}.
}
See \bold{References} and the references therein.

Standard errors are calculated using the sandwich (Godambe) information matrix.
}

\value{
	Returns a list where \code{par} are the estimated parameters, \code{LL} is the value of the maximized
	log-likelihood, \code{TIC} is the Takeuchi Information Criterion and \code{SE} are the standard errors.
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
# The following examples provide the fitting
# results of the air quality data recorded in
# the city center of Leeds, UK, analysed in
# Beranger and Padoan (2015).
################################################

## Load datsets
data(pollution)


## Dataset PM10-NO-SO2 (PNS)

\dontrun{
alik(PNS,model="Pairwise",c(1,1,1,1),trace=2,sig=2)
alik(PNS,model="Husler",rep(1,3),trace=2,sig=2)
alik(PNS,model="Dirichlet",rep(0.1,3),trace=2,sig=2)
alik(PNS,model="Extremalt",c(-0.5,-0.4,-0.5,1),c=0.01,trace=2,sig=2)
alik(PNS,model="Asymmetric",c(rep(1.1,4),rep(0.1,9)),c=0.01,trace=2,sig=2)
}

## Dataset NO2-SO2-NO (NSN)

\dontrun{
alik(NSN,model="Pairwise",c(1,1,1,1),trace=2,sig=2)
alik(NSN,model="Husler",rep(1,3),trace=2,sig=2)
alik(NSN,model="Dirichlet",rep(0.1,3),trace=2,sig=2)
alik(NSN,model="Extremalt",c(-0.5,-0.4,-0.5,1),c=0.01,trace=2,sig=2)
alik(NSN,model="Asymmetric",c(rep(1.1,4),rep(0.1,9)),c=0.01,trace=2,sig=2)
}

## Dataset PM10-NO-NO2 (PNN)

\dontrun{
alik(PNN,model="Pairwise",c(1,1,1,1),trace=2,sig=2)
alik(PNN,model="Husler",rep(1,3),trace=2,sig=2)
alik(PNN,model="Dirichlet",rep(0.1,3),trace=2,sig=2)
alik(PNN,model="Extremalt",c(-0.5,-0.4,-0.5,1),c=0.01,trace=2,sig=2)
alik(PNN,model="Asymmetric",c(rep(1.1,4),rep(0.1,9)),c=0.01,trace=2,sig=2)
}

## Dataset PM10-NO-NO2-SO2 (PNNS)

\dontrun{
alik(PNNS,model="Pairwise",rep(1,choose(ncol(PNNS),2)+1),trace=2,sig=2)
alik(PNNS,model="Husler",rep(1,choose(ncol(PNNS),2)),trace=2,sig=2)
alik(PNNS,model="Dirichlet",rep(1,ncol(PNNS)),trace=2,sig=2)
}

}

\keyword{optimize}
