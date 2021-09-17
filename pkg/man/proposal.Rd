\name{proposal}
\alias{proposal}

\title{Proposal distribution for parametric models}

\description{Density of the proposal distribution \code{q(cur.par, prop.par)} and
	random generator of the MCMC algorithm for parametric models. }

\usage{
	proposal(model, type = c("r", "d"), cur.par, prop.par, MCpar, log = TRUE)
}

\arguments{
	\item{model}{The parametric model considered. Values can be \code{model="Pairwise"}, \code{model="Husler"},
		\code{model="Dirichlet"}, \code{model="Extremalt"} and \code{model="Asymmetric"} respectively for the
		Pairwise Beta, Hulser-Reiss, Tilted Dirichlet, Extremal-t and Asymmetric Logistic.}
	\item{type}{One of the character strings "r" or "d" representing random generation and density for the
		Asymmetric Logistic model.}
	\item{cur.par}{Vector representing the current state of the chain. See \bold{Details}.}
	\item{prop.par}{Vector representing the candidate parameters. See \bold{Details}.}
	\item{MCpar}{A list made of a single element: MC MC parameter. Re-centering parameters for the
		proposal distribution.}
	\item{log}{Logical; Only used if \code{type=="d"} in order to obtain the log-density.
		\code{TRUE} is the default.}
}

\details{
	\itemize{
		\item For the \code{Pairwise Beta model}, \code{cur.par} and \code{prop.par} are of size \code{choose(dim,2)+1}.
	The components \code{prop.par[i]} of the proposal distribution are generated independently from
	the log-normal distribution.
	\code{prop.par = rlnorm(length(cur.par), meanlog=log(cur.par), sdlog = rep(MCpar$sdlog, length(cur.par)))};
		\item For the \code{Husler-Reiss model}, \code{cur.par} and \code{prop.par} are of size \code{choose(dim,2)}.
	The components \code{prop.par[i]} of the proposal distribution are generated independently from
	the log-normal distribution.
	\code{prop.par = rlnorm(length(cur.par), meanlog=log(cur.par), sdlog = rep(MCpar$sdlog, length(cur.par)))};
		\item For the \bold{Tilted Dirichlet model}, \code{cur.par} and \code{prop.par} are of size \code{dim}.
	The components \code{prop.par[i]} of the proposal distribution are generated independently from
	the log-normal distribution.
	\code{prop.par = rlnorm(length(cur.par), meanlog=log(cur.par), sdlog = rep(MCpar$sdlog, length(cur.par)))};
		\item For the \bold{Extremal-t model}, \code{cur.par} and \code{prop.par} are of size \code{choose(dim,2)+1}.
	The components \code{prop.par[i]} of the proposal distribution are generated independently from
	the square root of an inverse logit transformation of the normal distribution for the correlation
	parameters \code{rho} and from the log transformation of the normal distribution for the global dependence
	parameter (the degreee of freedom \code{mu});
		\item For the \bold{Asymmetric Logistic model}, \code{cur.par} and \code{prop.par} are
	of size \eqn{2^{d-1}(d+2)-(2d-1)}{2^{d-1}(d+2)-(2d-1)}.
	The components \code{prop.par[i]} of the proposal distribution are generated independently from
	the log transformation of the normal distribution (minus 1) for the parameters \code{alpha}
	and from the inverse logit transformation of the parameters \code{beta};
	}
}

\value{
	Either the (log)-density of the proposal \code{prop.par}, given \code{cur.par} (if \code{type=="d"})
	or a proposal parameter (a vector), if \code{type=="r".}
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://mypage.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com/};
}
\examples{

proposal("Pairwise", "r", rep(1,4), MCpar=0.35)
proposal("Pairwise", "d", rep(1,4), rep(1.2,4), MCpar=0.35)

proposal("Husler", "r", rep(1,4), MCpar=0.35)
proposal("Husler", "d", rep(1,4), c(1.2,4), MCpar=0.35)

proposal("Dirichlet", "r", rep(1,3), MCpar=0.35)
proposal("Dirichlet", "d", rep(1,3), c(1.2,3), MCpar=0.35)

proposal("Extremalt", "r", rep(0.5,3), MCpar=0.35)
proposal("Extremalt", "d", c(rep(0.9,3),3), c(rep(0.8,3),2), MCpar=0.35)

proposal("Asymmetric", "r", c(rep(1.1,4),rep(0.1,9)), MCpar=0.35)
proposal("Asymmetric", "d", c(rep(1.1,4),rep(0.1,9)), c(rep(1.2,4),rep(0.1,9)), 0.35)

}

\keyword{models}
