####################################################################################
### Authors: Boris Beranger and Simone Padoan        	 		                 ###
### 				   				                                             ###
###	Emails: borisberanger@gmail.com, simone.padoan@unibocconi.it				 ###
### 																			 ###
###	Institutions: Department of Decision Sciences, University Bocconi of Milan	 ###
### School of Mathematics and Statistics, University of New South Wales 		 ###
### 																			 ###
### File name: Densities.r	                 							         ###
### 																			 ###
### Description:                                  		              			 ###
### This file provides the Angular densities of extremal dependence models: 	 ###
### 1) Pairwise Beta (PB), Tilted Dirichlet (TD) and Husler-Reiss (HR)           ###
###     with only mass on the interior of the simplex)                           ###
### 2) Asymmetric logistic (AL), Extremal-t (ET) and extremal skew-t (EST)       ###
###     with mass on all subsets of the simplex)				                 ###
### It also provides the likelihood and log-likelihood functions				 ###
### 																			 ###
### Last change: 06/08/2019                         		  					 ###
### 																			 ###
####################################################################################

### Estimation and generation from angular density

angular <- function(data, model, n, dep, asy, alpha, beta, df, seed, k, nsim, plot=TRUE, nw=100){
  
  w <- seq(0.00001, .99999, length=nw)
  
  # Simulation of synthetic data
  if(missing(data)){
    if(missing(model) || missing(n) || missing(dep) || missing(seed) ){stop("model, n. dep and missing must be specified when data is not provided")}
    # warning("Data not provided and generated according to the parameters provided")
    
    if(any(model==c("log", "alog", "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"))){
      set.seed(seed)
      data <- rbvevd(n=n, dep = dep, asy, alpha, beta, model = model, mar1 = c(1,1,1))
      Atrue <- abvevd(w, dep=dep, asy, alpha, beta, model=model) # True pickands dependence function
      htrue <- hbvevd(1-w, dep=dep, asy, alpha, beta, model=model, half=TRUE) # True angular density
      if(any(model==c("alog","aneglog"))){Htrue <- (1-asy)/2}
      if(model=="amix"){Htrue <- c(1-alpha-beta, 1-alpha-2*beta)/2}
    }else if(model=="Extremalt"){
      set.seed(seed)
      data <- rExtDep(n=n, model=model, par=c(dep, df), angular=FALSE, mar=c(1,1,1), num=5e+5)
      Atrue <- rep(NA,nw)
      for(i in 1:nw){Atrue[i] <- index.ExtDep(object="pickands", model="ET", par=c(dep,df), x=c(w[i],1-w[i]))}
      htrue <- dExtDep(x=cbind(w,1-w), method="Parametric", model=model, par=c(dep,df), angular=TRUE, c=0, log=FALSE, vectorial=TRUE)/2
      Htrue <- dExtDep(x=cbind(c(1,0),c(0,1)), method="Parametric", model=model, par=c(dep,df), angular=TRUE, c=0, log=FALSE, vectorial=TRUE)/2
    }
  }else{
    if(!is.matrix(data) || ncol(data)!=2){stop("data must be a matrix with 2 columns")}
    n <- nrow(data)
    model <- dep <- seed <- NULL
    Atrue <- htrue <- 0
    warning("True Pickands function and angular density not computed as data is provided and true model is unsure")
  }
  
  if(missing(k)){stop("k, the polynomial order must be specified")}
  if(missing(nsim)){stop("nsim, the number of of generation from the estimated angular density must be specified")}  
  
  # Compute the angles:
  wdata <- data[,1]/rowSums(data)
  
  # Estimate the pickands function:
  Aest <- beed(data, cbind(w, 1-w), 2, 'md', 'emp', k=k, plot=FALSE)
  beta <- Aest$beta
  
  # Compute the angular density in the interior
  hest <- sapply(1:nw, function(i, w, beta) dh(w[i], beta), w, beta)
  
  # Compute the angular measure
  Hest <- sapply(1:nw, function(i, w, beta) ph(w[i], beta), w, beta)
  
  # Compute the point masses on the corners
  p0 <- (1+k*(beta[2]-1))/2
  p1 <- (1-k*(1-beta[k]))/2
  
  # simulate from the semiparametric angular density
  wsim <- rh(nsim, beta)
  
  if(plot){
    par(mai = c(0.85, 0.85, 0.1, 0.1), mgp = c(2.8, 1, 0), cex.axis=2, cex.lab=2)
    hist(wsim[wsim!=0 & wsim!=1], freq=FALSE, ylim=c(0,3.5), xlim=c(0,1), xlab='w', ylab='h(w)', main="",col="gray")
    lines(w, hest, col=1, lty=2, lwd=2.5)
    if(is.vector(htrue)){ lines(w[2:99], htrue[2:99],col=2, lty=1, lwd=1.5)}
    points(1,sum(wsim==1)/nsim, pch=19, cex=2)
    points(0,sum(wsim==0)/nsim, pch=19, cex=2)
    if(any(model==c('alog','aneglog','amix','Extremalt'))){
      points(0, Htrue[1] , pch=21, cex=2,col=2, lwd=2)
      points(1, Htrue[2], pch=21, cex=2,col=2, lwd=2)
    }
  }
  
  return(list(model=model, n=n, dep=dep, data=data, Aest=Aest, hest=hest, Hest=Hest, p0=p0, p1=p1, wsim=wsim, Atrue=Atrue, htrue=htrue))
  
}
