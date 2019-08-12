\name{prior}
\alias{prior}

\title{Prior parameter distribution for parametric models.}
\description{Random generation from the prior distribution for extremal parametric models or density evaluation of the extremal parametric models.}
\usage{
	prior(model, type = c("r", "d"), n, par, Hpar, log, dimData)
}

\arguments{
	\item{model}{The parametric model considered. Values can be \code{model="Pairwise"}, \code{model="Husler"},
		\code{model="Dirichlet"}, \code{model="Extremalt"} and \code{model="Asymmetric"} respectively for the
		Pairwise Beta, Hulser-Reiss, Tilted Dirichlet, Extremal-t and Asymmetric Logistic.}
	\item{type}{One of the character strings "r" or "d" representing random generation and density of the
		 model considered.}
	\item{n}{The number of parameters to be generated. Only used if \code{type=="r"}.}
	\item{par}{The values of the parameters.
		Only used if \code{type=="d"}. The first elements correspond to the parameters \code{alpha} and
		the last parameters are the \code{beta} parameters. See \bold{Details}.}
	\item{Hpar}{A list of hyper-parameters. See \bold{Details}.}
	\item{log}{Logical; Only used if \code{type=="d"} in order to obtain the log-density.
		\code{TRUE} is the default.}
	\item{dimData}{The dimension of the simplex.}
}

\details{
	\itemize{
		\item For the \bold{Pairwise Beta model}, the parameters components are independent, log-normal.
	The vector of parameters is of size \code{choose(dim,2)+1} with positive components. The first elements are the
	pairiwse dependence parameters \code{b} and the last one is the global dependence parameter \code{alpha}.
	The list of hyper-parameters should be of the form
	\code{mean.alpha=, mean.beta=, sd.alpha=, sd.beta=};
		\item For the \bold{Husler-Reiss model}, the parameters components are independent, log-normal.
	The vector of parameters is of size \code{choose(dim,2)+1} with positive components. The elements correspond
	to the \code{lambda} parameter. The list of hyper-parameters should be of the form \code{mean.lambda=, sd.lambda=};
		\item For the \bold{Dirichlet model}, the parameters' components are independent, log-normal.
	The vector of parameters is of size \code{dimData} with positive components. The elements correspond to the \code{alpha} parameter. The list of hyper-parameters should be of the form \code{mean.alpha=, sd.alpha=};
		\item For the \bold{Extremal-t model}, the parameters' components are independent, logit-squared for \code{rho} and log-normal for \code{mu}. The vector of parameters is of size \code{dimData} with positive components. The first elements correspond	to the correlation parameters \code{rho} and the last parameter is the global dependence parameter \code{mu}. The list of hyper-parameters should be of the form \code{mean.rho=, mean.mu=, sd.rho=, sd.mu=};
		\item For the \bold{Asymmetric Logistic model}, the parameters' components are independent, log(+1)-normal for \code{alpha} and logit for \code{beta}.
		The vector of parameters is of size
		\code{2^{dimData-1}(dimData+2)-(2 dimData+1)}{2^{dimData-1}(dimData+2)-(2 dimData+1)} with positive components.
		The list of hyper-parameters should be of the form \code{mean.alpha=, mean.beta=, sd.alpha=, sd.beta=}.
	}
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.unibocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{http://www.borisberanger.com};
}

\value{
	If \code{type=="r"}, a matrix with \eqn{n} rows containing a random parameter sample generated under the prior
	is returned, the (log)-density is returned if \code{type=="d"}.
}

\examples{
MCpar <- 0.35
Hpar.pb <- 	list(mean.alpha=0, mean.beta=3,sd.alpha=3, sd.beta=3)
Hpar.hr <- list(mean.lambda=0, sd.lambda=3)
Hpar.di <- list(mean.alpha=0, sd.alpha=3)
Hpar.et <- list(mean.rho=0, mean.mu=3,sd.rho=3, sd.mu=3)
Hpar.alm <- list(mean.alpha=0, mean.beta=0, sd.alpha=3, sd.beta=3)

prior(model="Pairwise", type="r", n=5, Hpar=Hpar.pb, dimData=3)
prior(model="Pairwise", type="d", par=rep(1,choose(4,2)+1), Hpar=Hpar.pb, log=TRUE, dimData=3)

prior(model="Husler", type="r", n=5, Hpar=Hpar.hr, dimData=3)
prior(model="Husler", type="d", par=rep(1,choose(4,2)), Hpar=Hpar.hr, log=TRUE, dimData=3)

prior(model="Dirichlet", type="r", n=5, Hpar=Hpar.di, dimData=3)
prior(model="Dirichlet", type="d", par=rep(1,3), Hpar=Hpar.di, log=TRUE, dimData=3)

prior(model="Extremalt", type="r", n=5, Hpar=Hpar.et, dimData=3)
prior(model="Extremalt", type="d", par=c(rep(0.1,3),4), Hpar=Hpar.et, log=TRUE, dimData=3)

prior(model="Asymmetric", type="r", n=5, Hpar=Hpar.alm, dimData=3)
prior(model="Asymmetric", type="d", par=c(rep(2,4),rep(0.7,9)), Hpar=Hpar.alm, log=TRUE, dimData=3)

}

\keyword{models}
