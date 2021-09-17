\name{posteriorMCMC}
\alias{posteriorMCMC}

\title{MCMC sampler for parametric spectral measures}
\description{Generates a sample from the posterior distribution for the parameters and computes the posterior mean,
	component-wise variance and BIC.}
\usage{
posteriorMCMC(Nsim, Nbin=0, Hpar, MCpar, dat, par.start=NULL,
              show.progress=floor(seq(1,Nsim, length.out=20)),
	            seed=NULL, kind="Mersenne-Twister", save=FALSE,
              name.save=NULL, save.directory = "~",
              name.dat="", model, c=NULL)
}

\arguments{
	\item{Nsim}{Total number of iterations to perform.}
	\item{Nbin}{Length of the burn-in period.}
	\item{Hpar}{A vector of hyper-parameters. See \link{prior}.}
	\item{MCpar}{MC MC parameter.See \link{proposal}.}
	\item{dat}{Angular dataset. Each row corresponds to coordinates in the simplex.}
	\item{par.start}{Starting point for the MC MC sample.}
	\item{show.progress}{A vector of integers containing the times (iteration numbers) at which a message showing
		progression will be printed on the standard output.}
	\item{seed}{The seed to be set via the routine \bold{set.seed}, see help of \bold{R} for details.}
	\item{kind}{The kind of random number generator. Default is \code{"Mersenne-Twister"}. See \bold{set.seed}
		for details.}
	\item{save}{Logical; if \code{save=TRUE} then the result is saved}
	\item{name.save}{A character string giving the name under which the result is to be saved. If \code{NULL} (default),
		nothing is saved. Otherwise the result is saved in file
		\code{paste(save.directory,"/", name.save,".rda",sep="")}.
		A "log" is also saved, named \code{paste(name.save, ".log", sep="")}, in file \code{paste(save.directory,"/",
		name.log, ".rda", sep="")}.}
	\item{save.directory}{A character string giving the directory where the result is to be saved (without trailing slash).}
	\item{name.dat}{A character string naming the dataset used for inference. Default is "".}
	\item{model}{A character string. Possible models are \code{"Pairiwse"}, \code{"Husler"}, \code{"Dirichlet"},
		\code{"Extremalt"} or \code{"Asymmetric"}. See \code{details}.}
	\item{c}{A real value in \eqn{[0,1]}{[0,1]}, providing the decision rule to allocate a data point to a subset of the simplex.
		Only required for the Extremal-t and Asymmetric Logistic models.}
}

\details{
	When \code{model="Pairwise"} the Pairiwse Beta model is selected and \code{prior.pb}, \code{proposal.pb}, \code{pb.Hpar}
	, \code{pb.MCpar} are considered. Similarly \code{model="Husler"} selects the Husler-Reiss model,
	\code{model="Dirichlet"} the Tilted Dirichlet model, \code{model="Extremalt"} the Extremal-t and \code{model="Asymmetric"}
	the Asymmetric Logistic model and the functions associated to these models.
}

\value{
	A list made of
	\item{stored.vales}{A \eqn{(Nsim-Nbin)*d}{(Nsim-Nbin)*d} matrix, where \eqn{d}{d} is the dimension of the parameter space}
	\item{llh}{A vector of size \eqn{(Nsim-Nbin)}{(Nsim-Nbin)} containing the log-likelihoods evaluadted at each parameter
		of the posterior sample.}
	\item{lprior}{A vector of size \eqn{(Nsim-Nbin)}{(Nsim-Nbin)} containing the logarithm of the prior densities evaluated
		at each parameter of the posterior sample.}
	\item{elapsed}{The time elapsed, as given by \code{porc.time} between the start and end of the run.}
	\item{Nsim}{The same as the passed argument.}
	\item{Nbin}{Idem.}
	\item{n.accept}{The total number of accepted proposals.}
	\item{n.accept.kept}{The number of accepted proposals after the burn-in period. }
	\item{emp.mean}{The estimated posterior parameters mean.}
	\item{emp.sd}{The empirical posterior sample standard deviation.}
	\item{BIC}{The Bayesian Information Criteria.}
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://mypage.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com/};
}


\examples{

################################################
# The following examples provide the results of
# the approximate bayesian analysis in Table 1.1
# of the paper Beranger and Padoan (2015)
################################################

## Load datasets :
data(pollution)
Nsim <- 50e+4
Nbin <- 30e+4
MCpar <- 0.35
Hpar.pb <- 	list(mean.alpha=0, mean.beta=3,sd.alpha=3, sd.beta=3)
Hpar.hr <- list(mean.lambda=0, sd.lambda=3)
Hpar.di <- list(mean.alpha=0, sd.alpha=3)
Hpar.et <- list(mean.rho=0, mean.mu=3,sd.rho=3, sd.mu=3)

## Using the PNS dataset
if (interactive()){
est.pb.PNS <- posteriorMCMC(Nsim, Nbin, Hpar.pb, MCpar, PNS, seed=14342, model='Pairwise')
est.pb.PNS$emp.mean
est.pb.PNS$emp.sd
est.pb.PNS$BIC

est.hr.PNS <- posteriorMCMC(Nsim, Nbin, Hpar.hr, MCpar, PNS, seed=14342, model='Husler')
est.hr.PNS$emp.mean
est.hr.PNS$emp.sd
est.hr.PNS$BIC

est.di.PNS <- posteriorMCMC(Nsim, Nbin, Hpar.di, MCpar, PNS, seed=14342, model='Dirichlet')
est.di.PNS$emp.mean
est.di.PNS$emp.sd
est.di.PNS$BIC

est.et.PNS <- posteriorMCMC(Nsim, Nbin, Hpar.et, MCpar, PNS, seed=14342, model='Extremalt',c=0.1)
est.et.PNS$emp.mean
est.et.PNS$emp.sd
est.et.PNS$BIC
}

## Using the NSN dataset

if (interactive()){
est.pb.NSN <- posteriorMCMC(Nsim, Nbin, Hpar.pb, MCpar, NSN, seed=14342, model='Pairwise')
est.pb.NSN$emp.mean
est.pb.NSN$emp.sd
est.pb.NSN$BIC

est.hr.NSN <- posteriorMCMC(Nsim, Nbin, Hpar.hr, MCpar, NSN, seed=14342, model='Husler')
est.hr.NSN$emp.mean
est.hr.NSN$emp.sd
est.hr.NSN$BIC

est.di.NSN <- posteriorMCMC(Nsim, Nbin, Hpar.di, MCpar, NSN, seed=14342, model='Dirichlet')
est.di.NSN$emp.mean
est.di.NSN$emp.sd
est.di.NSN$BIC

est.et.NSN <- posteriorMCMC(Nsim, Nbin, Hpar.et, MCpar, NSN, seed=14342, model='Extremalt',c=0.1)
est.et.NSN$emp.mean
est.et.NSN$emp.sd
est.et.NSN$BIC
}

## Using the PNN dataset

if (interactive()){

est.pb.PNN <- posteriorMCMC(Nsim, Nbin, Hpar.pb, MCpar, PNN, seed=14342, model='Pairwise')
est.pb.PNN$emp.mean
est.pb.PNN$emp.sd
est.pb.PNN$BIC

est.hr.PNN <- posteriorMCMC(Nsim, Nbin, Hpar.hr, MCpar, PNN, seed=14342, model='Husler')
est.hr.PNN$emp.mean
est.hr.PNN$emp.sd
est.hr.PNN$BIC

est.di.PNN <- posteriorMCMC(Nsim, Nbin, Hpar.di, MCpar, PNN, seed=14342, model='Dirichlet')
est.di.PNN$emp.mean
est.di.PNN$emp.sd
est.di.PNN$BIC

est.et.PNN <- posteriorMCMC(Nsim, Nbin, Hpar.et, MCpar, PNN, seed=14342, model='Extremalt',c=0.1)
est.et.PNN$emp.mean
est.et.PNN$emp.sd
est.et.PNN$BIC
}

################################################
# The following examples provide the results of
# the approximate bayesian analysis in Table 1.2
# of the paper Beranger and Padoan (2015)
################################################


# Using the PNNS dataset

if (interactive()){

est.pb.PNNS <- posteriorMCMC(Nsim, Nbin, Hpar.pb, MCpar, PNNS, seed=14342, model='Pairwise')
est.pb.PNNS$BIC

est.hr.PNNS <- posteriorMCMC(Nsim, Nbin, Hpar.hr, MCpar, PNNS, seed=14342, model='Husler')
est.hr.PNNS$BIC

est.di.PNNS <- posteriorMCMC(Nsim, Nbin, Hpar.di, MCpar, PNNS, seed=14342, model='Dirichlet')
est.di.PNNS$BIC
}

}

\keyword{models}
