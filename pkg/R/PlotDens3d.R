################################################################################################
### Authors: Boris Beranger and Simone Padoan        	 									 ###
### 																							 ###	
###	Emails: borisberanger@gmail.com, simone.padoan@unibocconi.it								 ###
### 																							 ###
###	Institutions: Department of Decision Sciences, University Bocconi of Milan				 ###
### School of Mathematics and Statistics, University of New South Wales 						 ###
### 																							 ###
### File name: PlotDens3d.r				      							             	     ###
### 																							 ###
### Description:                                  							      		     ###
### This file provides plots of the trivariate angular densities				 			 	 ###
### for the Pairwise Beta, Dirichlet, H??sler-Reiss, Asym Logisitic and Extremal-t models. 	 ###
### 																							 ###
### Last change: 31/08/2016                         		  									 ###
### 																							 ###
################################################################################################



AngDensPlot <- function(model='Pairwise', para=c(2,4,15,1), log=TRUE, data=NULL, contour=TRUE, 
labels=c(expression(w[1]),expression(w[3]),expression(w[2])),cex.dat=1, cex.lab=1, cex.cont=1){

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
			lf[i] = dens(x=x[i,], model=model, par=para, c=0, log=log, vectorial=TRUE)			
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
			quant[i] = dens(x=data[i,], model=model, par=para, c=0, log=log, vectorial=TRUE)
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