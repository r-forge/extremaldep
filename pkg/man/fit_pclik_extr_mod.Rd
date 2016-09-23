\name{fit_pclik_extr_mod}
\alias{fit_pclik_extr_mod}

\title{Fit extremal dependence models using pairwise composite likelihood}
\description{Estimate the parameters of the Husler-Reiss, Extremal-$t$ and Extremal Skew-$t$ models using pairwise composite likelihood, for up to \eqn{4}{4} dimensional datasets. }

\usage{
	fit_pclik_extr_mod(model, data, parastart, trace)
}

\arguments{
	\item{model}{A string with the name of the model: \code{"hr"}, \code{"Extremalt"} or \code{"Skewt"}. }
	\item{data}{A \code{data.frame} or \code{matrix} obejct with up to \eqn{4}{4} columns.}
	\item{parastart}{A vector containing the initial parameter values. See \bold{Details}. }	
	\item{trace}{A non-negative integer. If positive, tracing information on the progress of the optimization is produced.
		See \bold{?optim} for details.}
}

\details{
	
	Data must be marginally on unit Frechet scale.
	
	If \code{model="hr"} then the vector of initial values is made of \code{choose(d,2)} positive parameters, \code{d=2,3}.
	If \code{model="Extremalt"} then the vector of initial values is made of \code{choose(d,2)} dependence parameters and a degree of freedom, \code{d=2,3}.
	If \code{model="Skewt"} then the vector of initial values is made of \code{choose(d,2)} dependence parameters, \code{d} shape (or skewness) parameters and a degree of freedom, \code{d=2,3}.
	
	In the case of bivariate data the regular likelihood estimation is performed.	
	
}

\value{
	Returns the vector of estimated parameters and the value of the pairwise composite log-likelihood.
}

\references{
	
	Beranger, B. and Padoan, S. A. (2015).
	Extreme dependence models, chapater of the book \emph{Extreme Value Modeling and Risk Analysis: Methods and Applications},
	\bold{Chapman Hall/CRC}.
	
	Beranger, B., Padoan, S. A. and Sisson, S. A. (2016).
	Models for extremal dependence derived from skew-symmetric families.
	\emph{Scandinavian Journal of Statistics}, \bold{To appear}.	
	
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.bocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
}

\examples{

## Reproduce the real data analysis from
## Beranger et al. (2016), Section 5.

data(Wind)

## Vector of starting values
p0 <- c(rep(0.5,3),1)

### CLOU CLAY SALL

\dontrun{		
ext1 <- fit_pclik_extr_mod("Extremalt", CLOU.CLAY.SALL, p0, 2)
est.ext1 <- round(ext1$par,2)
p01 <- c(ext1$par[1:3],rep(0,3),ext1$par[4])
skewt1 <- fit_pclik_extr_mod("Skewt", CLOU.CLAY.SALL, p01, 2)
est.skewt1 <- round(skewt1$par,2)
}

### CLOU CLAY PAUL

\dontrun{		
ext2 <- fit_pclik_extr_mod("Extremalt", CLOU.CLAY.PAUL, p0, 2)
est.ext2 <- round(ext2$par,2)
p02 <- c(ext2$par[1:3],rep(0,3),ext2$par[4])
skewt2 <- fit_pclik_extr_mod("Skewt", CLOU.CLAY.PAUL, p02, 2)
est.skewt2 <- round(skewt2$par,2)
}

### CLAY SALL PAUL

\dontrun{	
ext3 <- fit_pclik_extr_mod("Extremalt", CLAY.SALL.PAUL, p0, 2)
est.ext3 <- round(ext3$par,2)
p03 <- c(ext3$par[1:3],rep(0,3),ext3$par[4])
skewt3 <- fit_pclik_extr_mod("Skewt", CLAY.SALL.PAUL, p03, 2)
est.skewt3 <- round(skewt3$par,2)
}

### CLAY SALL PAUL

\dontrun{	
ext4 <- fit_pclik_extr_mod("Extremalt", CLOU.SALL.PAUL, p0, 2)
est.ext4 <- round(ext4$par,2)
p04 <- c(ext4$par[1:3],rep(0,3),ext4$par[4])
skewt4 <- fit_pclik_extr_mod("Skewt", CLOU.SALL.PAUL, p04, 2)
est.skewt4 <- round(skewt4$par,2)
}

}




