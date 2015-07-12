################################################################################################
### Authors: Boris Beranger and Simone Padoan        	 									 ###
### 																						 ###	
###	Emails: borisberanger@gmail.com, simone.padoan@unibocconi.it							 ###
### 																						 ###
###	Institutions: Department of Decision Sciences, University Bocconi of Milan				 ###
### School of Mathematics and Statistics, University of New South Wales 					 ###
### 																						 ###
### File name: Exponent.r   	                 							             	 ###
### 																						 ###
### Description:                                  							      		     ###
### This file provides the Exponent function and calculates the probability of excedding     ###
### some threshold for the extremal dependence models.										 ###
### So far only the Husler-Reiss model is inlcuded.										     ###
### 																						 ###
### Last change: 11/07/2015                         		  								 ###
### 																						 ###
################################################################################################

###
### Husler-Reiss model
###

# Exponent function (see Beranger and Padoan (2014), p.16)
# Dimensions d=2,3,4,5

exponent.hr <- function(X, par, dist=FALSE){

	expo.hr.2d <- function(X,par){
		x=X[1]; y=X[2];
		lambda12 <- par;
		pt1 <- 1/x * pnorm(lambda12+log(y/x)/(2*lambda12))
		pt2 <- 1/y * pnorm(lambda12+log(x/y)/(2*lambda12))
		return( pt1 + pt2 )
	}

	expo.hr.3d <- function(X,par){
		x=X[1]; y=X[2]; z=X[3];
		lambda12 <- par[1]; lambda13 <- par[2]; lambda23 <- par[3];
		S1 <- matrix(c(4*lambda12^2, 2*(lambda12^2 + lambda13^2 - lambda23^2) , 2*(lambda12^2 + lambda13^2 - lambda23^2) , 4*lambda13^2),nrow=2)
		S2 <- matrix(c(4*lambda12^2, 2*(lambda12^2 + lambda23^2 - lambda13^2) , 2*(lambda12^2 + lambda23^2 - lambda13^2) , 4*lambda23^2),nrow=2)
		S3 <- matrix(c(4*lambda13^2, 2*(lambda13^2 + lambda23^2 - lambda12^2) , 2*(lambda13^2 + lambda23^2 - lambda12^2) , 4*lambda23^2),nrow=2)
		pt1 <- 1/x * pmvnorm(lower=rep(-Inf,2),upper=c(log(y/x)+2*lambda12^2, log(z/x)+2*lambda13^2),sigma = S1 )
		pt2 <- 1/y * pmvnorm(lower=rep(-Inf,2),upper=c(log(x/y)+2*lambda12^2, log(z/y)+2*lambda23^2),sigma = S2 )
		pt3 <- 1/z * pmvnorm(lower=rep(-Inf,2),upper=c(log(x/z)+2*lambda13^2, log(y/z)+2*lambda23^2),sigma = S3 )
		return(pt1 + pt2 + pt3 )
	}

	expo.hr.4d <- function(X,par){
		x=X[1]; y=X[2]; z=X[3]; w=X[4];
		lambda12 <- par[1]; lambda13 <- par[2]; lambda14 <- par[3]; lambda23 <- par[4]; lambda24 <- par[5]; lambda34 <- par[6];

		R1 <- matrix(c(4*lambda12^2, 2*(lambda12^2 + lambda13^2 -lambda23^2), 2*(lambda12^2 + lambda14^2 -lambda24^2),
		2*(lambda12^2 + lambda13^2 -lambda23^2), 4*lambda13^2, 2*(lambda13^2 + lambda14^2 -lambda34^2),
		2*(lambda12^2 + lambda14^2 -lambda24^2), 2*(lambda13^2 + lambda14^2 - lambda34^2), 4*lambda14^2),nrow=3)
	
		R2 <- matrix(c(4*lambda12^2, 2*(lambda12^2 + lambda23^2 -lambda13^2), 2*(lambda12^2 + lambda24^2 -lambda14^2),
		2*(lambda12^2 + lambda23^2 -lambda13^2), 4*lambda23^2, 2*(lambda23^2 + lambda24^2 -lambda34^2),
		2*(lambda12^2 + lambda24^2 -lambda14^2), 2*(lambda23^2 + lambda24^2 - lambda34^2), 4*lambda24^2),nrow=3)
	
		R3 <- matrix(c(4*lambda13^2, 2*(lambda13^2 + lambda23^2 -lambda12^2), 2*(lambda13^2 + lambda34^2 -lambda14^2),
		2*(lambda13^2 + lambda23^2 -lambda12^2), 4*lambda23^2, 2*(lambda23^2 + lambda34^2 -lambda24^2),
		2*(lambda13^2 + lambda34^2 -lambda14^2), 2*(lambda23^2 + lambda34^2 -lambda24^2), 4*lambda34^2),nrow=3)

		R4 <- matrix(c(4*lambda14^2, 2*(lambda14^2 + lambda24^2 -lambda12^2), 2*(lambda14^2 + lambda34^2 -lambda13^2),
		2*(lambda14^2 + lambda24^2 -lambda12^2), 4*lambda24^2, 2*(lambda24^2 + lambda34^2 -lambda23^2),
		2*(lambda14^2 + lambda34^2 -lambda13^2), 2*(lambda24^2 + lambda34^2 -lambda23^2), 4*lambda34^2),nrow=3)

		pt1 <- 1/x *pmvnorm(lower=rep(-Inf,3), upper=c(2*lambda12^2+log(y/x),2*lambda13^2+log(z/x),2*lambda14^2+log(w/x) ), sigma = R1)
		pt2 <- 1/y *pmvnorm(lower=rep(-Inf,3), upper=c(2*lambda12^2+log(x/y),2*lambda23^2+log(z/y),2*lambda24^2+log(w/y) ), sigma = R2)
		pt3 <- 1/z *pmvnorm(lower=rep(-Inf,3), upper=c(2*lambda13^2+log(x/z),2*lambda23^2+log(y/z),2*lambda34^2+log(w/z) ), sigma = R3)
		pt4 <- 1/w *pmvnorm(lower=rep(-Inf,3), upper=c(2*lambda14^2+log(x/w),2*lambda24^2+log(y/w),2*lambda34^2+log(z/w) ), sigma = R4)	
		return(pt1 + pt2 + pt3 + pt4)
	}

	expo.hr.5d <- function(X,par){
		x=X[1]; y=X[2]; z=X[3]; s=X[4]; t=X[5]
		lambda12 <- par[1]; lambda13 <- par[2]; lambda14 <- par[3]; lambda15 <- par[4]; lambda23 <- par[5];
		lambda24 <- par[6]; lambda25 <- par[7]; lambda34 <- par[8]; lambda35 <- par[9]; lambda45 <- par[10];
		R1 <- matrix(c(4*lambda12^2, 2*(lambda12^2+lambda13^2-lambda23^2), 2*(lambda12^2+lambda14^2-lambda24^2), 2*(lambda12^2+lambda15^2-lambda25^2),
		2*(lambda12^2+lambda13^2-lambda23^2), 4*lambda13^2, 2*(lambda13^2+lambda14^2-lambda34^2), 2*(lambda13^2+lambda15^2-lambda35^2),
		2*(lambda12^2+lambda14^2-lambda24^2), 2*(lambda13^2+lambda14^2-lambda34^2), 4*lambda14^2, 2*(lambda14^2+lambda15^2-lambda45^2),
		2*(lambda12^2+lambda15^2-lambda25^2), 2*(lambda13^2+lambda15^2-lambda35^2), 2*(lambda14^2+lambda15^2-lambda45^2), 4*lambda15^2
		),ncol=4)
	
		R2 <- matrix(c(4*lambda12^2, 2*(lambda12^2+lambda23^2-lambda13^2), 2*(lambda12^2+lambda24^2-lambda14^2), 2*(lambda12^2+lambda25^2-lambda15^2),
		2*(lambda12^2+lambda23^2-lambda13^2), 4*lambda23^2, 2*(lambda23^2+lambda24^2-lambda34^2), 2*(lambda23^2+lambda25^2-lambda35^2),
		2*(lambda12^2+lambda24^2-lambda14^2), 2*(lambda23^2+lambda24^2-lambda34^2), 4*lambda24^2, 2*(lambda24^2+lambda25^2-lambda45^2),
		2*(lambda12^2+lambda25^2-lambda15^2), 2*(lambda23^2+lambda25^2-lambda35^2), 2*(lambda24^2+lambda25^2-lambda45^2), 4*lambda25^2
		),ncol=4)
	
		R3 <- matrix(c(4*lambda13^2, 2*(lambda13^2+lambda23^2-lambda12^2), 2*(lambda13^2+lambda34^2-lambda14^2), 2*(lambda13^2+lambda35^2-lambda15^2),
		2*(lambda13^2+lambda23^2-lambda12^2), 4*lambda23^2, 2*(lambda23^2+lambda34^2-lambda24^2), 2*(lambda23^2+lambda35^2-lambda25^2),
		2*(lambda13^2+lambda34^2-lambda14^2), 2*(lambda23^2+lambda34^2-lambda24^2), 4*lambda34^2, 2*(lambda34^2+lambda35^2-lambda45^2),
		2*(lambda13^2+lambda35^2-lambda15^2), 2*(lambda23^2+lambda35^2-lambda25^2), 2*(lambda34^2+lambda35^2-lambda45^2), 4*lambda35^2
		),ncol=4)

		R4 <- matrix(c(4*lambda14^2, 2*(lambda14^2+lambda24^2-lambda12^2), 2*(lambda14^2+lambda34^2-lambda13^2), 2*(lambda14^2+lambda45^2-lambda15^2),
		2*(lambda14^2+lambda24^2-lambda12^2), 4*lambda24^2, 2*(lambda24^2+lambda34^2-lambda23^2), 2*(lambda24^2+lambda45^2-lambda25^2),
		2*(lambda14^2+lambda34^2-lambda13^2), 2*(lambda24^2+lambda34^2-lambda23^2), 4*lambda34^2, 2*(lambda34^2+lambda45^2-lambda35^2),
		2*(lambda14^2+lambda45^2-lambda15^2), 2*(lambda24^2+lambda45^2-lambda25^2), 2*(lambda34^2+lambda45^2-lambda35^2), 4*lambda45^2
		),ncol=4)

		R5 <- matrix(c(4*lambda15^2, 2*(lambda15^2+lambda25^2-lambda12^2), 2*(lambda15^2+lambda35^2-lambda13^2), 2*(lambda15^2+lambda45^2-lambda14^2),
		2*(lambda15^2+lambda24^2-lambda12^2), 4*lambda25^2, 2*(lambda25^2+lambda35^2-lambda23^2), 2*(lambda25^2+lambda45^2-lambda24^2),
		2*(lambda15^2+lambda34^2-lambda13^2), 2*(lambda25^2+lambda35^2-lambda23^2), 4*lambda35^2, 2*(lambda35^2+lambda45^2-lambda45^2),
		2*(lambda15^2+lambda45^2-lambda14^2), 2*(lambda25^2+lambda45^2-lambda24^2), 2*(lambda35^2+lambda45^2-lambda34^2), 4*lambda45^2
		),ncol=4)

		pt1 <- 1/x *pmvnorm(lower=rep(-Inf,4), upper=c(2*lambda12^2+log(y/x),2*lambda13^2+log(z/x),2*lambda14^2+log(s/x),2*lambda15^2+log(t/x) ), sigma = R1)
		pt2 <- 1/y *pmvnorm(lower=rep(-Inf,4), upper=c(2*lambda12^2+log(x/y),2*lambda23^2+log(z/y),2*lambda24^2+log(s/y),2*lambda25^2+log(t/y) ), sigma = R2)
		pt3 <- 1/z *pmvnorm(lower=rep(-Inf,4), upper=c(2*lambda13^2+log(x/z),2*lambda23^2+log(y/z),2*lambda34^2+log(s/z),2*lambda35^2+log(t/z) ), sigma = R3)
		pt4 <- 1/s *pmvnorm(lower=rep(-Inf,4), upper=c(2*lambda14^2+log(x/s),2*lambda24^2+log(y/s),2*lambda34^2+log(z/s),2*lambda45^2+log(t/s) ), sigma = R4)
		pt5 <- 1/t *pmvnorm(lower=rep(-Inf,4), upper=c(2*lambda15^2+log(x/t),2*lambda25^2+log(y/t),2*lambda35^2+log(z/t),2*lambda45^2+log(s/t) ), sigma = R5)
		return(pt1 + pt2 + pt3 + pt4 + pt5)
	}

	dim=length(X)
	npar=length(par)
	if(any(par<=0)){break('wrong value of parameters')}	
	if(dim>5){break('exponent and distribution functions not available for dimensions greater than 5')}
	if(dim==2 && npar==1){res <- expo.hr.2d(X,par)} # in 2d there is only 1 parameter
	if(dim==3 && npar==3){res <- expo.hr.3d(X,par)} # in 3d there are 3 parameters
	if(dim==4 && npar==6){res <- expo.hr.4d(X,par)} # in 4d there are 6 parameters
	if(dim==5 && npar==10){res <- expo.hr.5d(X,par)} # in 5d there are 10 parameters
	
	if(dist==TRUE){return(exp(-res))}else{
		return(res)
	}	
}


excess.dist.hr <- function(X, par){

	excess.dist.hr.2d <- function(X,par){
		x=X[1]; y=X[2];
		lambda12 <- par;
		pt1 <- 1/x * pnorm(-lambda12-log(y/x)/(2*lambda12))
		pt2 <- 1/y * pnorm(-lambda12-log(x/y)/(2*lambda12))
		return( pt1 + pt2 )
	}

	excess.dist.hr.3d <- function(X,par){
		x=X[1]; y=X[2]; z=X[3];
		lambda12 <- par[1]; lambda13 <- par[2]; lambda23 <- par[3];
		S1 <- matrix(c(4*lambda12^2, 2*(lambda12^2 + lambda13^2 - lambda23^2) , 2*(lambda12^2 + lambda13^2 - lambda23^2) , 4*lambda13^2),nrow=2)
		S2 <- matrix(c(4*lambda12^2, 2*(lambda12^2 + lambda23^2 - lambda13^2) , 2*(lambda12^2 + lambda23^2 - lambda13^2) , 4*lambda23^2),nrow=2)
		S3 <- matrix(c(4*lambda13^2, 2*(lambda13^2 + lambda23^2 - lambda12^2) , 2*(lambda13^2 + lambda23^2 - lambda12^2) , 4*lambda23^2),nrow=2)
		pt1 <- 1/x * pmvnorm(lower=rep(-Inf,2),upper=c(-log(y/x)-2*lambda12^2, -log(z/x)-2*lambda13^2),sigma = S1 )
		pt2 <- 1/y * pmvnorm(lower=rep(-Inf,2),upper=c(-log(x/y)-2*lambda12^2, -log(z/y)-2*lambda23^2),sigma = S2 )
		pt3 <- 1/z * pmvnorm(lower=rep(-Inf,2),upper=c(-log(x/z)-2*lambda13^2, -log(y/z)-2*lambda23^2),sigma = S3 )
		return(pt1 + pt2 + pt3 )
	}

	excess.dist.hr.4d <- function(X,par){
		x=X[1]; y=X[2]; z=X[3]; w=X[4];
		lambda12 <- par[1]; lambda13 <- par[2]; lambda14 <- par[3]; lambda23 <- par[4]; lambda24 <- par[5]; lambda34 <- par[6];
		R1 <- matrix(c(4*lambda12^2, 2*(lambda12^2 + lambda13^2 -lambda23^2), 2*(lambda12^2 + lambda14^2 -lambda24^2),
		2*(lambda12^2 + lambda13^2 -lambda23^2), 4*lambda13^2, 2*(lambda13^2 + lambda14^2 -lambda34^2),
		2*(lambda12^2 + lambda14^2 -lambda24^2), 2*(lambda13^2 + lambda14^2 - lambda34^2), 4*lambda14^2),nrow=3)
	
		R2 <- matrix(c(4*lambda12^2, 2*(lambda12^2 + lambda23^2 -lambda13^2), 2*(lambda12^2 + lambda24^2 -lambda14^2),
		2*(lambda12^2 + lambda23^2 -lambda13^2), 4*lambda23^2, 2*(lambda23^2 + lambda24^2 -lambda34^2),
		2*(lambda12^2 + lambda24^2 -lambda14^2), 2*(lambda23^2 + lambda24^2 - lambda34^2), 4*lambda24^2),nrow=3)
	
		R3 <- matrix(c(4*lambda13^2, 2*(lambda13^2 + lambda23^2 -lambda12^2), 2*(lambda13^2 + lambda34^2 -lambda14^2),
		2*(lambda13^2 + lambda23^2 -lambda12^2), 4*lambda23^2, 2*(lambda23^2 + lambda34^2 -lambda24^2),
		2*(lambda13^2 + lambda34^2 -lambda14^2), 2*(lambda23^2 + lambda34^2 -lambda24^2), 4*lambda34^2),nrow=3)

		R4 <- matrix(c(4*lambda14^2, 2*(lambda14^2 + lambda24^2 -lambda12^2), 2*(lambda14^2 + lambda34^2 -lambda13^2),
		2*(lambda14^2 + lambda24^2 -lambda12^2), 4*lambda24^2, 2*(lambda24^2 + lambda34^2 -lambda23^2),
		2*(lambda14^2 + lambda34^2 -lambda13^2), 2*(lambda24^2 + lambda34^2 -lambda23^2), 4*lambda34^2),nrow=3)

		pt1 <- 1/x *pmvnorm(lower=rep(-Inf,3), upper=c(-2*lambda12^2-log(y/x),-2*lambda13^2-log(z/x),-2*lambda14^2-log(w/x) ), sigma = R1)
		pt2 <- 1/y *pmvnorm(lower=rep(-Inf,3), upper=c(-2*lambda12^2-log(x/y),-2*lambda23^2-log(z/y),-2*lambda24^2-log(w/y) ), sigma = R2)
		pt3 <- 1/z *pmvnorm(lower=rep(-Inf,3), upper=c(-2*lambda13^2-log(x/z),-2*lambda23^2-log(y/z),-2*lambda34^2-log(w/z) ), sigma = R3)
		pt4 <- 1/w *pmvnorm(lower=rep(-Inf,3), upper=c(-2*lambda14^2-log(x/w),-2*lambda24^2-log(y/w),-2*lambda34^2-log(z/w) ), sigma = R4)
		return(pt1 + pt2 + pt3 + pt4)
	}

	excess.dist.hr.5d <- function(X,par){
		x=X[1]; y=X[2]; z=X[3]; s=X[4]; t=X[5]
		lambda12 <- par[1]; lambda13 <- par[2]; lambda14 <- par[3]; lambda15 <- par[4]; lambda23 <- par[5];
		lambda24 <- par[6]; lambda25 <- par[7]; lambda34 <- par[8]; lambda35 <- par[9]; lambda45 <- par[10];
		R1 <- matrix(c(4*lambda12^2, 2*(lambda12^2+lambda13^2-lambda23^2), 2*(lambda12^2+lambda14^2-lambda24^2), 2*(lambda12^2+lambda15^2-lambda25^2),
		2*(lambda12^2+lambda13^2-lambda23^2), 4*lambda13^2, 2*(lambda13^2+lambda14^2-lambda34^2), 2*(lambda13^2+lambda15^2-lambda35^2),
		2*(lambda12^2+lambda14^2-lambda24^2), 2*(lambda13^2+lambda14^2-lambda34^2), 4*lambda14^2, 2*(lambda14^2+lambda15^2-lambda45^2),
		2*(lambda12^2+lambda15^2-lambda25^2), 2*(lambda13^2+lambda15^2-lambda35^2), 2*(lambda14^2+lambda15^2-lambda45^2), 4*lambda15^2
		),ncol=4)
	
		R2 <- matrix(c(4*lambda12^2, 2*(lambda12^2+lambda23^2-lambda13^2), 2*(lambda12^2+lambda24^2-lambda14^2), 2*(lambda12^2+lambda25^2-lambda15^2),
		2*(lambda12^2+lambda23^2-lambda13^2), 4*lambda23^2, 2*(lambda23^2+lambda24^2-lambda34^2), 2*(lambda23^2+lambda25^2-lambda35^2),
		2*(lambda12^2+lambda24^2-lambda14^2), 2*(lambda23^2+lambda24^2-lambda34^2), 4*lambda24^2, 2*(lambda24^2+lambda25^2-lambda45^2),
		2*(lambda12^2+lambda25^2-lambda15^2), 2*(lambda23^2+lambda25^2-lambda35^2), 2*(lambda24^2+lambda25^2-lambda45^2), 4*lambda25^2
		),ncol=4)
	
		R3 <- matrix(c(4*lambda13^2, 2*(lambda13^2+lambda23^2-lambda12^2), 2*(lambda13^2+lambda34^2-lambda14^2), 2*(lambda13^2+lambda35^2-lambda15^2),
		2*(lambda13^2+lambda23^2-lambda12^2), 4*lambda23^2, 2*(lambda23^2+lambda34^2-lambda24^2), 2*(lambda23^2+lambda35^2-lambda25^2),
		2*(lambda13^2+lambda34^2-lambda14^2), 2*(lambda23^2+lambda34^2-lambda24^2), 4*lambda34^2, 2*(lambda34^2+lambda35^2-lambda45^2),
		2*(lambda13^2+lambda35^2-lambda15^2), 2*(lambda23^2+lambda35^2-lambda25^2), 2*(lambda34^2+lambda35^2-lambda45^2), 4*lambda35^2
		),ncol=4)

		R4 <- matrix(c(4*lambda14^2, 2*(lambda14^2+lambda24^2-lambda12^2), 2*(lambda14^2+lambda34^2-lambda13^2), 2*(lambda14^2+lambda45^2-lambda15^2),
		2*(lambda14^2+lambda24^2-lambda12^2), 4*lambda24^2, 2*(lambda24^2+lambda34^2-lambda23^2), 2*(lambda24^2+lambda45^2-lambda25^2),
		2*(lambda14^2+lambda34^2-lambda13^2), 2*(lambda24^2+lambda34^2-lambda23^2), 4*lambda34^2, 2*(lambda34^2+lambda45^2-lambda35^2),
		2*(lambda14^2+lambda45^2-lambda15^2), 2*(lambda24^2+lambda45^2-lambda25^2), 2*(lambda34^2+lambda45^2-lambda35^2), 4*lambda45^2
		),ncol=4)

		R5 <- matrix(c(4*lambda15^2, 2*(lambda15^2+lambda25^2-lambda12^2), 2*(lambda15^2+lambda35^2-lambda13^2), 2*(lambda15^2+lambda45^2-lambda14^2),
		2*(lambda15^2+lambda24^2-lambda12^2), 4*lambda25^2, 2*(lambda25^2+lambda35^2-lambda23^2), 2*(lambda25^2+lambda45^2-lambda24^2),
		2*(lambda15^2+lambda34^2-lambda13^2), 2*(lambda25^2+lambda35^2-lambda23^2), 4*lambda35^2, 2*(lambda35^2+lambda45^2-lambda45^2),
		2*(lambda15^2+lambda45^2-lambda14^2), 2*(lambda25^2+lambda45^2-lambda24^2), 2*(lambda35^2+lambda45^2-lambda34^2), 4*lambda45^2
		),ncol=4)

		pt1 <- 1/x *pmvnorm(lower=rep(-Inf,4), upper=c(-2*lambda12^2+log(y/x),-2*lambda13^2+log(z/x),-2*lambda14^2+log(s/x),-2*lambda15^2+log(t/x) ), sigma = R1)
		pt2 <- 1/y *pmvnorm(lower=rep(-Inf,4), upper=c(-2*lambda12^2+log(x/y),-2*lambda23^2+log(z/y),-2*lambda24^2+log(s/y),-2*lambda25^2+log(t/y) ), sigma = R2)
		pt3 <- 1/z *pmvnorm(lower=rep(-Inf,4), upper=c(-2*lambda13^2+log(x/z),-2*lambda23^2+log(y/z),-2*lambda34^2+log(s/z),-2*lambda35^2+log(t/z) ), sigma = R3)
		pt4 <- 1/s *pmvnorm(lower=rep(-Inf,4), upper=c(-2*lambda14^2+log(x/s),-2*lambda24^2+log(y/s),-2*lambda34^2+log(z/s),-2*lambda45^2+log(t/s) ), sigma = R4)
		pt5 <- 1/t *pmvnorm(lower=rep(-Inf,4), upper=c(-2*lambda15^2+log(x/t),-2*lambda25^2+log(y/t),-2*lambda35^2+log(z/t),-2*lambda45^2+log(s/t) ), sigma = R5)
		return(pt1 + pt2 + pt3 + pt4 + pt5)

	}

	dim=length(X)
	npar=length(par)
	if(any(par<=0)){break('wrong value of parameters')}
	if(dim>5){break('exponent function not available for dimensions greater than 5')}
	if(dim==2 && npar==1){return(excess.dist.hr.2d(X,par))} # in 2d there is only 1 parameter
	if(dim==3 && npar==3){return(excess.dist.hr.3d(X,par))} # in 3d there are 3 parameters
	if(dim==4 && npar==6){return(excess.dist.hr.4d(X,par))} # in 4d there are 6 parameters
	if(dim==5 && npar==10){return(excess.dist.hr.5d(X,par))} # in 5d there are 10 parameters
}





