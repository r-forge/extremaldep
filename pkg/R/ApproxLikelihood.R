################################################################################################
### Authors: Boris Beranger and Simone Padoan        	 									 ###
### 																							 ###	
###	Emails: borisberanger@gmail.com, simone.padoan@unibocconi.it								 ###
### 																							 ###
###	Institutions: Department of Decision Sciences, University Bocconi of Milan				 ###
### School of Mathematics and Statistics, University of New South Wales 						 ###
### 																							 ###
### File name: ApproxLikelihood.r              							             	     ###
### 																							 ###
### Description:                                  							      		     ###
### This file provides the Approximate Likelihood estimation 				 				 ###
### for the Pairwise Beta, Dirichlet, H??sler-Reiss, Asym Logisitic and Extremal-t models. 	 ###
### Also calculates the standard errors using the Godambe information matrix.				 ###
### 																							 ###
### Last change: 10/12/2014                         		  									 ###
### 																							 ###
################################################################################################

###############################################################################
## Trace of Matrix (identical to matrixcalc package for example) 
###############################################################################

matrix.trace <- function (x) 
{
    if (ncol(x)!=nrow(x)) 
        stop("argument x is not a square matrix")
    return(sum(diag(x)))
}

###############################################################################
## Matrix square root - taken from Stephen Lake 
## http://www5.biostat.wustl.edu/s-news/s-news-archive/200109/msg00067.html
###############################################################################

matrix.sqrt <- function(A)
{
	if (length(A)==1)
    return(sqrt(A))
	sva <- svd(A)
	if (min(sva$d)>=0)
    Asqrt <- sva$u %*% diag(sqrt(sva$d)) %*% t(sva$v)
	else
    stop("matrix square root is not defined")
	return(Asqrt)
}

alik <- function(data, model, parastart, c=NULL, trace=0, sig=3){

	optimfun<-function(para){
		return(dens(x=data, model=model, par=para, c=c, log=TRUE, vectorial=FALSE))
	}
	est.para=optim(parastart, optimfun,method="BFGS",control=list(maxit=500, trace=trace, fnscale=-1)) 
	
	LogLik <- est.para$value # log likelihood
	param.est <- est.para$par # Parameters estimates

	n <- nrow(data)
	s <- 0
	score <- matrix(nrow=n,ncol=length(param.est))
	for(i in 1:n){
		stepJ <- function(para){
			dens(x=data[i,], model=model, par=para, c=c, log=TRUE, vectorial=FALSE)
		}
		score[i,] = jacobian(stepJ,param.est)
		s = s + hessian(stepJ,param.est)	

	} 
	K=var(score); # variability matrix
	J=-s; # sensitivity matrix
	TIC = 2*matrix.trace(K %*% solve(J))-2*LogLik; # TIC
	SE = diag(matrix.sqrt(solve((J %*% solve(K) %*% J)/n ))) # Standard errors
		
	return(list(par=round(param.est,sig), LL=round(LogLik,sig), TIC=round(TIC,sig), SE=round(SE,sig) ))
}