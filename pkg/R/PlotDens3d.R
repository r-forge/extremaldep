################################################################################################
### Authors: Boris Beranger and Simone Padoan        	 									 ###
### 																						 ###	
###	Emails: borisberanger@gmail.com, simone.padoan@unibocconi.it							 ###
### 																						 ###
###	Institutions: Department of Decision Sciences, University Bocconi of Milan				 ###
### School of Mathematics and Statistics, University of New South Wales 					 ###
### 																						 ###
### File name: PlotDens3d.r				      							             	     ###
### 																						 ###
### Description:                                  							      		     ###
### This file provides plots of the trivariate angular densities				 		 	 ###
### for the Pairwise Beta, Dirichlet, H??sler-Reiss, Asym Logisitic and Extremal-t models. 	 ###
### 																						 ###
### Last change: 11/07/2015                         		  								 ###
### 																						 ###
################################################################################################



AngDensPlot <- function(model='Pairwise', para=c(2,4,15,1), log=TRUE, data=NULL, contour=TRUE, 
labels=c(expression(w[1]),expression(w[3]),expression(w[2])),cex.dat=1, cex.lab=1, cex.cont=1){

### SPECIFIC DENSITIES FROM AL AND ET MODEL (WE ONLY PLOT THE MASS ON THE INTERIOR OF THE SIMPLEX)

	interior_alm_3d_2 <- function(w,alpha,beta){
	
		k <- c(1,2,3)
	
		part1 <- (alpha-1) * (2*alpha-1)
		part2 <- prod( beta[k]^alpha * w[k]^(-alpha-1) )
		part3 <- sum( (beta[k]/w[k])^alpha )^(1/alpha-3)
	
		return(part1*part2*part3)
	}
  
	interior_et_d <- function(w,rho,mu){ # mass on the interior of the d-dim simplex
	  p<-length(w);
	  k=p-1;
	  w.tilde<-rep(0,k);
	  Sigma<-diag(k);
	  
	  for(i in 1:k){
	    w.tilde[i] <- ((w[i+1]/w[1])^(1/mu)-rho[i])*sqrt(mu+1);
	    for(j in i:k){
	      if(i==j){Sigma[i,j]=Sigma[j,i]=1-rho[i]^2}else{
	        Sigma[i,j]=Sigma[j,i]=(rho[sum(k:(k-i+1))+j-i]-rho[i]*rho[j])
	      }
	    }
	  }
	  
	  if(sum(eigen(Sigma)$values<0)>=1){return(1e-50)} #Check if matrix is positive definite
	  
	  deriv = (w[-1]/w[1])^(1/mu-1)/mu*sqrt(mu+1)
	  return(dmvt(w.tilde,rep(0,k),sigma=Sigma,df=mu+1,log=FALSE)*w[1]^(-p-1)*prod(deriv) )	
	  
	}

	dens_al_2 <- function(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), alpha=c(1.2),beta=rep(0.3,3), log=FALSE, vectorial=TRUE){
		xvect = as.double(as.vector(t(x)))
   		if (is.vector(x)) {
        	n = as.integer(1)
       		if(round(sum(x),7) !=1){ stop("Data is not angular") }
        	result <- interior_alm_3d_2(x,alpha,beta)
    	}
    	else {
       		n = as.integer(nrow(x))
    		if (sum(apply(x,1,sum)) != n){ stop("Data is not angular") }
	    	if (vectorial) {
	       		result = double(n)
	        	result <- apply(x,1,function(y){interior_alm_3d_2(y,alpha,beta)})
	    	} else { # vectorial=FALSE mean we return the likelihood function
    	    	result = as.double(1)
        		result <- prod(apply(x,1,function(y){interior_alm_3d_2(y,alpha,beta)})) 
    		}
    	}	
    	if(log)
    		return(log(result))
    	else return(result)
	}

	dens_et_2 <- function(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), rho=rep(0.1,3),mu=2, log=FALSE, vectorial=TRUE){
		xvect = as.double(as.vector(t(x)))
	    if (is.vector(x)) {
 	       dim = as.integer(length(x))
 	       n = as.integer(1)
 	       if(round(sum(x),7) != 1){ stop("Data is not angular")}
 	   	result <- interior_et_d(x,rho,mu)
 	   }
 	   else {
 	       dim = as.integer(ncol(x))
 	       n = as.integer(nrow(x))
    	
			if (sum(apply(x,1,sum)) != n){ stop("Data is not angular") }
    
 	   	if (vectorial) {
 	   	    result = double(n)
 	   	    result <- apply(x,1,function(y){interior_et_d(y,rho,mu)}) 
 	   	} else { # vectorial=FALSE mean we return the likelihood function
 	   	    result = as.double(1)
 	   	    result <- prod(apply(x,1,function(y){interior_et_d(y,rho,mu)}))
 	  		}
 	  	}	
 	   if(log)
 	   	return(log(result))
 	   else return(result)
	}

	RectTri<-function(x,model,para,log){
	  if(is.vector(x)){x<-matrix(x,nrow=1)}
	  n<-nrow(x)
	  d<-ncol(x)
	  ind<-(rowSums(x>0)==d)*(rowSums(x)<1) # Give which coordinates have both components greater than zero and the sum of the components is less than 1
	  ind [ind == 0] <- NA
	  x<-cbind(x,1-rowSums(x))
 	 lf<-vector("numeric")
  
 	 for(i in 1:nrow(x)){
 	 	if(prod(x[i,]>=0)==1){
			if(model=='Extremalt'){lf[i] = dens_et_2(x=x[i,], rho=para[1:3], mu=para[4], log=log, vectorial=TRUE)}
 	 		if(model=='Asymmetric'){lf[i] = dens_al_2(x=x[i,], alpha=para[1], beta=para[2:4], log=log, vectorial=TRUE)}
			if(model!='Extremalt' && model!='Asymmetric'){lf[i] = dens(x=x[i,], model=model, par=para, log=log, vectorial=TRUE)}
 	 	} else {lf[i]=0}
 	 }
  
 	 f<-lf*ind    
 	 return(f)
	}

	EquiTri <-function(x,model,para,log){
	  y <- x
	  y[,1] <- x[,1]- x[,2]/sqrt(3)
	  y[,2] <- x[,2]*2/sqrt(3)
	  fx <- RectTri(x=y, model=model, para=para, log=log)/(sqrt(3)/2) ## adjust with Jacobian
	  equil.ind <- !((sqrt(3)*x[,1] - x[,2] <= 0) & x[,1]<=1/2 | (sqrt(3)*x[,1] + x[,2] >= sqrt(3)) & x[,1]>=1/2)  ## indicator for equilateral triangle
	  equil.ind [equil.ind == 0] <- NA
	  return(fx*equil.ind)
	
	}

	x1 <- seq(0,1, length=301)
	x2 <- seq(0,1, length=301)
	xx <- as.matrix(expand.grid(x1,x2))

	equi <- EquiTri(x=xx, model= model, para=para, log=log)
	dec <- seq (from=0.1, to=0.9, by=0.1)

	if(!is.null(data)){
		quant=0
		for(i in 1:nrow(data)){
			if(model=='Extremalt'){ quant[i] = dens_et_2(x=data[i,], rho=para[1:3], mu=para[4], log=log, vectorial=TRUE)}
		  	if(model=='Asymmetric'){ quant[i] = dens_al_2(x=data[i,], alpha=para[1], beta=para[2:4], log=log, vectorial=TRUE)}	
			if(model!='Extremalt' && model!='Asymmetric'){quant[i] = dens(x=data[i,], model=model, par=para, log=log, vectorial=TRUE)}
		}
		deciles <- c(min(equi,na.rm=TRUE),quantile(quant,dec,na.rm=TRUE)/(sqrt(3)/2),max(equi,na.rm=TRUE))
	}else{
		deciles <- c(min(equi,na.rm=TRUE),quantile(equi,dec,na.rm=TRUE),max(equi,na.rm=TRUE))
	}	

	image(x1,x2, matrix(equi, nrow=length(x1), ncol=length(x2)), asp=1, breaks=deciles, col=rev(heat.colors(length(deciles)-1)),
        axes=FALSE,xlab="",ylab="",xlim=c(-0.05,1.06),ylim=c(-0.08,0.99))
	if(!is.null(data)){ points(data[,1] + 0.5*data[,2],sqrt(3)/2*data[,2],pch=16,cex=cex.dat)}
	text(c(0.94,0.5,0.06),c(-0.05,0.905,-0.05),labels=labels,cex=cex.lab)
	segments(c(0,0,0.5),c(0,0,sqrt(3)/2),c(1,0.5,1),c(0,sqrt(3)/2,0))
	if(contour==TRUE){ contour(x1,x2, matrix(equi, nrow=length(x1), ncol=length(x2)),levels=deciles, 
                             labels=round(deciles,3), labcex=cex.cont, nlevels=15,add=TRUE)}
}