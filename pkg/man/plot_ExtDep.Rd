\name{plot_ExtDep}
\alias{plot_ExtDep}

\title{Graphical summaries of parametric representations of extremal dependence.}
\description{This function displays the angular density, Pickands dependence function and return levels for bivariate and trivariate extreme values models.}

\usage{
plot_ExtDep(object="angular", model, par, log=TRUE, data=NULL, contour=TRUE, 
            style, labels, cex.dat=1, cex.lab=1, cex.cont=1, 
            Q.fix, Q.range, Q.range0, cond=FALSE, ...)
}

\arguments{
	\item{object}{A character string indicating which graphical summary to plot. Takes value \code{"angular"} (default) \code{"pickands"} or\code{"returns"}.}
	\item{model}{A string with the name of the model considered. Takes value \code{"PB"} (Pairwise Beta), \code{"HR"} (Husler-Reiss), \code{"ET"} (Extremal-t), \code{"EST"} (Extremal Skew-t), \code{TD} (Tilted Dirichlet) or \code{AL} (Asymmetric Logistic) when evaluating the angular density. Restricted to \code{"HR"}, \code{"ET"} and \code{"EST"} for the Pickands dependence function.}
	\item{par}{A vector representing the parameters of the model.}
	\item{log}{A logical value specifying if the log density is computed. Required when \code{object="angular"}.}
	\item{data}{A matrix representing angular data to be added to the density plot. Required when \code{object="angular"}.}
	\item{contour}{A logical value; if \code{TRUE} the contour levels are displayed. Required for trivariate models only.}
	\item{style}{A character string indicating the plotting style of the data. Takes value \code{"hist"} or \code{"ticks"}. See \bold{details}. }
  \item{labels}{A vector of character strings indicating the labels. Must be of length \eqn{1}{1} for bivariate models and \eqn{3}{3} for trivariate models.}
  \item{cex.dat}{A positive real indicating the size of the data points. Required for the trivariate angular density.}
  \item{cex.lab}{A positive real indicating the size of the labels.}
  \item{cex.cont}{A positive real indicating the size of the contour labels.}
  \item{Q.fix}{A vector of length the dimension of the model, indicating some fixed quantiles to compute joint return levels. Must contain \code{NA} (maximum 2) for components whose quantiles are allowed to vary. Required when \code{object="returns"}.}
  \item{Q.range}{A vector or matrix indicating quantile values on the unit Frechet scale, for the components that are allowed to vary. Must be a vector or a one-column matrix if there is one \code{NA} in \code{Q.fix}. Must be a two-column matrix if there are two \code{NA}s in \code{Q.fix}. Required when \code{object="returns"}.}
  \item{Q.range0}{A object of the same format as \code{Q.range}, corresponding to the quantiles on the original scale. Required when \code{object="returns"}.}
  \item{cond}{A logical value; if \code{TRUE} then conditional return levels are computed where the conditional event is the fixed event. Default if \code{FALSE}. Required when \code{object="returns"}.}
  \item{...}{Additional graphical arguments for the \code{hist()} and \code{plot()} functions respectively used to compute the bivariate angular density and Pickands dependence function as well as for the \code{plot()} and \code{contour()} functionswhen \code{object="returns"}.}
}

\details{
  The angular density is computed using the function \code{\link{dExtDep}} with arguments \code{method="Parametric"} and \code{angular=TRUE}. The Pickands dependence function is computed using the function \code{\link{index.ExtDep}} with argument \code{object="pickands"}.
  
  When displaying the bivariate angular density and some data are provided (a 2-column matrix is specified for \code{data}), there is the choice to summarise the data using a histogram (\code{style="hist"}) or to display the observations using tick marks (\code{style="ticks"}). 
  
  When displaying return levels, there are two possibilities: univariate and bivariate return levels. Since the model dimensions are restricted to a maximum of three, in that case, aunivariate return level corresponds to fixing two components while a bivariate return level fixes only one component. The choice of the fixed component is decided by the position of the \code{NA} value(s) in the \code{Q.fix} argument. If \code{par} is a vector then the corresponding return level(s) are printed. However if \code{par} is a matrix then the return level(s) are evaluated for each value of the parameter vector and the mean, and empirical \eqn{95\%}{95\%} empirical interval are displayed. Typically this is used when posterior samples are available. When \code{par} is a matrix with only two rows, resulting plots may not provide much information.  
  
  When contours are displayed, levels are chosen to be the deciles. 
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{https://faculty.unibocconi.it/simonepadoan/};
	Boris Beranger, \email{borisberanger@gmail.com}
	\url{https://www.borisberanger.com};
}


\seealso{\code{\link{dExtDep}}, \code{\link{index.ExtDep}}. }

\examples{

data(pollution)

###############################
### Trivariate Husler-Reiss ###
###############################

if(interactive()){
  f.hr <- fExtDep(method="PPP", data=PNS, model="HR", par.start=rep(1,3))

  plot_ExtDep(object="angular", model="HR", par=f.hr$par, data=PNS, 
              labels=c(expression(PM[10]), expression(NO), expression(SO[2])), 
              cex.lab=2)
      
  plot_ExtDep(object="pickands", model="HR", par=f.hr$par, data=PNS, 
              labels=c(expression(PM[10]), expression(NO), expression(SO[2])), 
              cex.lab=2) # Takes time!      
}

###############################
### Bivariate Husler-Reiss ###
###############################            
            
PN <- na.omit(Leeds.frechet[,1:2])
PN <- cbind(PN, rowSums(PN))
PN <- PN[order(PN[,3], decreasing = TRUE)[1:100],]
PN <- PN[,1:2]/PN[,3]

f.hr2 <- fExtDep(method="PPP", data=PN, model = "HR", par.start = 1)
plot_ExtDep(model="HR", par=f.hr2$par, log=FALSE, data=PN, style="hist")
plot_ExtDep(model="HR", par=f.hr2$par, log=FALSE, data=PN, style="ticks") 
plot_ExtDep(object="pickands", model="HR", par=f.hr2$par) 

}

\keyword{models}
