################################################################################################
### Authors: Boris Beranger and Simone Padoan        	 									 ###
### 																							 ###	
###	Emails: borisberanger@gmail.com, simone.padoan@unibocconi.it								 ###
### 																							 ###
###	Institutions: Department of Decision Sciences, University Bocconi of Milan				 ###
### School of Mathematics and Statistics, University of New South Wales 						 ###
### 																							 ###
### File name: ApproxBayesian.r	                 							             	 ###
### 																							 ###
### Description:                                  							      		     ###
### This file provides the Approximate Bayesian procedure for the extremal dependence models ###
### Pairwise Beta, Dirichlet, Husler-Reiss, Asymmetric logistic and Extremal-t			     ###
### 																							 ###
### Last change: 15/12/2014                         		  									 ###
### 																							 ###
################################################################################################


prior <- function(model, type = c("r", "d"), n, par, Hpar, log, dimData){

	logit <- function(p){
		if(any(p<=0) || any(p>=1)){stop('p must be in [0,1]')}
		return( log(p) - log(1-p) )
	}

	invlogit <- function(x){
		return( 1/ ( 1 + exp(-x) ) )
	}

	###
	### Pairwise Beta model
	###

	# Small change compare to the function given in BMAmevt package: in type="r" ncol=3 changed into ncol = lengthPar-1

	prior.pb <- function (type = c("r", "d"), n, par, Hpar, log, dimData) {
	   	if (type == "r") {
	        p <- dimData
	        lengthPar <- choose(p, 2) + 1
	        alpha <- exp(rnorm(n, mean = Hpar$mean.alpha, sd = Hpar$sd.alpha))
	        beta <- exp(matrix(rnorm((lengthPar - 1) * n, mean = Hpar$mean.beta, sd = Hpar$sd.beta), ncol = lengthPar-1))
	        res <- cbind(alpha, beta)
	        return(res)
	    }
	    if (type == "d") {
	        lpar <- log(par)
	        ld1 <- dnorm(lpar[1], mean = Hpar$mean.alpha, sd = Hpar$sd.alpha, 
	            log = TRUE)
	        ld2 <- dnorm(lpar[-1], mean = Hpar$mean.beta, sd = Hpar$sd.beta, 
	            log = TRUE)
	        if (log) 
	            return(ld1 + sum(ld2))
	        else return(exp(ld1 + sum(ld2)))
	    }
	    stop("wrong 'type' argument")
	}

	###
	### Husler-Reiss model
	###


	prior.hr <- function(type = c("r", "d"), n, par, Hpar, log, dimData){
	    if (type == "r") {
	        p <- dimData
	        lengthPar <- choose(p, 2)
	        lambda <- exp(matrix(rnorm(lengthPar*n, mean = Hpar$mean.lambda, sd = Hpar$sd.lambda),ncol=lengthPar))
	        return(lambda)
	    }
	    if (type == "d") {
	        lpar <- log(par)
	        ld1 <- dnorm(lpar, mean = Hpar$mean.lambda, sd = Hpar$sd.lambda, log = TRUE)
	        if (log) 
	            return( sum(ld1) )
	        else return( exp( sum(ld1) ) )
	    }
	    stop("wrong 'type' argument")
	}

	###
	### Dirichlet model
	###


	prior.di <- function(type = c("r", "d"), n, par, Hpar, log, dimData){
	    if (type == "r") {
	        p <- dimData
	        lengthPar <- p
	        alpha <- exp(matrix(rnorm(lengthPar*n, mean = Hpar$mean.alpha, sd = Hpar$sd.alpha),ncol=lengthPar))
	        return(alpha)
	    }
	    if (type == "d") {
	        lpar <- log(par)
	        ld1 <- dnorm(lpar, mean = Hpar$mean.alpha, sd = Hpar$sd.alpha, log = TRUE)
	        if (log) 
	            return( sum(ld1) )
	        else return( exp( sum(ld1) ) )
	    }
	    stop("wrong 'type' argument")
	}

	###
	### Extremal-t
	###


	prior.et <- function(type = c("r", "d"), n, par, Hpar, log, dimData){
	    if (type == "r") {
	        p <- dimData
	        lengthPar <- choose(p, 2)+1
	        rho <- sqrt( invlogit( matrix(rnorm((lengthPar-1)*n, mean = Hpar$mean.rho, sd = Hpar$sd.rho),ncol=lengthPar-1)))
	        mu <- exp(rnorm(n, mean = Hpar$mean.mu, sd = Hpar$sd.mu))
	        res <- cbind(rho,mu)
	        return(res)
	    }
	    if (type == "d") {
	        lpar <- c( logit(par[-length(par)]^2), log(par[length(par)]) )
	        ld1 <- dnorm(lpar[-length(lpar)], mean = Hpar$mean.rho, sd = Hpar$sd.rho, log = TRUE )
	        ld2 <- dnorm(lpar[length(lpar)], mean = Hpar$mean.mu, sd = Hpar$sd.mu, log = TRUE )
	        if (log) 
	            return( sum(ld1) + ld2 )
	        else return( exp( sum(ld1) + ld2 ) )
	    }
	    stop("wrong 'type' argument")
	}

	###
	### Asymmetric Logistic
	###


	prior.al <- function(type = c("r", "d"), n, par, Hpar, log, dimData){
	    p <- dimData
		if (type == "r") {
	        if(p==2){
		        lengthPar <- 3
		        alpha <- exp(rnorm(n, mean = Hpar$mean.alpha, sd = Hpar$sd.alpha))+1
			    beta <- invlogit( matrix(rnorm(2*n, mean = Hpar$mean.beta, sd = Hpar$sd.beta),ncol=2))
		        res <- cbind(alpha,beta)
		        return(res)
			}	        
	        if(p==3){
		        lengthPar <- 13
		        alpha <- exp(matrix(rnorm(4*n, mean = Hpar$mean.alpha, sd = Hpar$sd.alpha),ncol=4))+1
			    beta <- invlogit( matrix(rnorm(9*n, mean = Hpar$mean.beta, sd = Hpar$sd.beta),ncol=9))
		        res <- cbind(alpha,beta)
		        return(res)
			}
	        if(p==4){
		        lengthPar <- 39
		        alpha <- exp(matrix(rnorm(11*n, mean = Hpar$mean.alpha, sd = Hpar$sd.alpha),ncol=11))+1
			    beta <- invlogit( matrix(rnorm(28*n, mean = Hpar$mean.beta, sd = Hpar$sd.beta),ncol=28))
		        res <- cbind(alpha,beta)
		        return(res)
			}
	    }
    	if (type == "d") {
    		if(p==2){
	    	    lpar <- c( log( par[1]-1), logit(par[2:3]) )
    		    ld1 <- dnorm(lpar[1], mean = Hpar$mean.alpha, sd = Hpar$sd.alpha, log = TRUE)
        		ld2 <- dnorm(lpar[2:3], mean = Hpar$mean.beta, sd = Hpar$sd.beta, log = TRUE)
        	}	
    		if(p==3){
	    	    lpar <- c( log( par[1:4]-1), logit(par[5:13]))
    		    ld1 <- dnorm(lpar[1:4], mean = Hpar$mean.alpha, sd = Hpar$sd.alpha, log = TRUE)
    	   	 	ld2 <- dnorm(lpar[5:13], mean = Hpar$mean.beta, sd = Hpar$sd.beta, log = TRUE)
    	    }
    		if(p==4){
		    lpar <- c( log( par[1:11]-1), logit(par[12:39]) )
    		    ld1 <- dnorm(lpar[1:11], mean = Hpar$mean.alpha, sd = Hpar$sd.alpha, log = TRUE)
   	 	    	ld2 <- dnorm(lpar[12:39], mean = Hpar$mean.beta, sd = Hpar$sd.beta, log = TRUE)
        	}        
        	if (log) 
        	    return( sum(ld1) + sum(ld2) )
        	else return( exp( sum(ld1) + sum(ld2) ) )
    	}
    	stop("wrong 'type' argument")
	}

	if(model=="Pairwise"){return( prior.pb(type, n, par, Hpar, log, dimData) )}
	if(model=="Husler"){return( prior.hr(type, n, par, Hpar, log, dimData) )}
	if(model=="Dirichlet"){return( prior.di(type, n, par, Hpar, log, dimData) )}
	if(model=="Extremalt"){return( prior.et(type, n, par, Hpar, log, dimData) )}
	if(model=="Asymmetric"){return( prior.al(type, n, par, Hpar, log, dimData) )}				

}


proposal <- function (model, type = c("r", "d"), cur.par, prop.par, MCpar, log = TRUE) { 

	logit <- function(p){
		if(any(p<0) || any(p>1)){stop('p must be in [0,1]')}
		return( log(p) - log(1-p) )
	}

	invlogit <- function(x){
		return( 1/ ( 1 + exp(-x) ) )
	}

	# Function from the BMAmevt package

	proposal.pb <- function (type = c("r", "d"), cur.par, prop.par, MCpar, log = TRUE) { 
	    sd <- rep(MCpar, length(cur.par))
	    transfo <- function(x) { log(x) }
	    invtransfo <- function(x) { exp(x) }
	    mean <- transfo(cur.par)
	    if (type == "r") {
	        return(invtransfo(rnorm(length(cur.par), mean = mean, sd = sd)))
	    }
	    if (type == "d") {
	        vect.res = sapply(1:length(prop.par), function(j) {
	            dnorm(transfo(prop.par[j]), mean = mean[j], sd = sd[j], log = TRUE)
	        })
	        return(ifelse(log, sum(vect.res), exp(sum(vect.res))))
	    }
	    stop("wrong type specification")
	}

	# Function from the BMAmevt package

	proposal.hr <- function (type = c("r", "d"), cur.par, prop.par, MCpar, log = TRUE) {
	    sd <- rep(MCpar, length(cur.par))
	    transfo <- function(x) { log(x) }
	    invtransfo <- function(x) { exp(x) }
	    mean <- transfo(cur.par)
	    if (type == "r") {
	        return(invtransfo(rnorm(length(cur.par), mean = mean, sd = sd)))        
	    }
	    if (type == "d") {
	        vect.res = sapply(1:length(prop.par), function(j) {
	            dnorm(transfo(prop.par[j]), mean = mean[j], sd = sd[j], log = TRUE)            
	        })
	        return(ifelse(log, sum(vect.res), exp(sum(vect.res))))
	    }
	    stop("wrong type specification")
	}

	proposal.di <- function (type = c("r", "d"), cur.par, prop.par, MCpar, log = TRUE) {
	    sd <- rep(MCpar, length(cur.par))
	    transfo <- function(x) {
	        log(x)
	    }
	    invtransfo <- function(x) {
	        exp(x)
	    }
	    mean <- transfo(cur.par)
	    if (type == "r") {
	        return(invtransfo(rnorm(length(cur.par), mean = mean, sd = sd)))        
	    }
	    if (type == "d") {
	        vect.res = sapply(1:length(prop.par), function(j) {
	            dnorm(transfo(prop.par[j]), mean = mean[j], sd = sd[j], log = TRUE)            
	        })
	        return(ifelse(log, sum(vect.res), exp(sum(vect.res))))
	    }
	    stop("wrong type specification")
	}

	proposal.et <- function (type = c("r", "d"), cur.par, prop.par, MCpar, log = TRUE) {
	    sd.rho <- rep(MCpar, length(cur.par)-1)
	    sd.df <-  MCpar
	    
	    transfo.rho <- function(x) { logit(x^2) }
	    invtransfo.rho <- function(x) { sign(x)*sqrt(invlogit(x)) }
	    transfo.df <- function(x) { log(x) }
	    invtransfo.df <- function(x) { exp(x) }

		mean.rho <- transfo.rho(cur.par[-length(cur.par)])
	    mean.df <- transfo.df(cur.par[length(cur.par)])
	    if (type == "r") {
	    	res <- c(invtransfo.rho(rnorm(length(cur.par)-1, mean = mean.rho, sd = sd.rho)), invtransfo.df(rnorm(1, mean = mean.df, sd = sd.df)))
	        return(res)        
	    }
	    if (type == "d") {
	        vect.res = c( sapply(1:(length(prop.par)-1), function(j) {dnorm(transfo.rho(prop.par[j]), mean = mean.rho[j], sd = sd.rho[j], log=TRUE)}), dnorm(transfo.df(prop.par[length(prop.par)]), mean = mean.df, sd = sd.df, log = TRUE))
	        return(ifelse(log, sum(vect.res), exp(sum(vect.res))))
	    }
	    stop("wrong type specification")
	}

	proposal.al <- function (type = c("r", "d"), cur.par, prop.par, MCpar, log = TRUE) {

	 	if(length(cur.par)==3){ 
		    sd.alpha <- rep(MCpar, 1)
		    sd.beta <-  rep(MCpar, 2)
	    }
	 	if(length(cur.par)==13){ 
		    sd.alpha <- rep(MCpar, 4)
		    sd.beta <-  rep(MCpar, 9)
	    }
	 	if(length(cur.par)==39){ 
		    sd.alpha <- rep(MCpar, 11)
		    sd.beta <-  rep(MCpar, 28)
	    }    

	    transfo.alpha <- function(x) { log(x-1) }
	    invtransfo.alpha <- function(x) { exp(x)+1 }    
	    transfo.beta <- function(x) { logit(x) }
	    invtransfo.beta <- function(x) { invlogit(x) }

		if(length(cur.par)==3){
			mean.alpha <- transfo.alpha(cur.par[1])
		    mean.beta <- transfo.beta(cur.par[2:3])
		}
		if(length(cur.par)==13){
			mean.alpha <- transfo.alpha(cur.par[1:4])
		    mean.beta <- transfo.beta(cur.par[5:13])
		}
		if(length(cur.par)==39){
			mean.alpha <- transfo.alpha(cur.par[1:11])
		    mean.beta <- transfo.beta(cur.par[12:39])
		}		
	    if (type == "r") {
	    	if(length(cur.par)==3){
		    	res <- c(invtransfo.alpha(rnorm(1, mean = mean.alpha, sd = sd.alpha)), invtransfo.beta(rnorm(2, mean = mean.beta, sd = sd.beta)))
		    }	
	    	if(length(cur.par)==13){
		    	res <- c(invtransfo.alpha(rnorm(4, mean = mean.alpha, sd = sd.alpha)), invtransfo.beta(rnorm(9, mean = mean.beta, sd = sd.beta)))
		    }
	    	if(length(cur.par)==39){
		    	res <- c(invtransfo.alpha(rnorm(11, mean = mean.alpha, sd = sd.alpha)), invtransfo.beta(rnorm(28, mean = mean.beta, sd = sd.beta)))
		    }
	        return(res)        
	    }
	    if (type == "d") {
	    	if(length(prop.par)==3){
		        vect.res = c( dnorm(transfo.alpha(prop.par[1]), mean = mean.alpha, sd = sd.alpha, log = TRUE), sapply(2:3, function(j) {dnorm(transfo.beta(prop.par[j]), mean = mean.beta[j-1], sd = sd.beta[j-1], log=TRUE)}) )
		    }
		    if(length(prop.par)==13){
		        vect.res = c( sapply(1:4, function(j) {dnorm(transfo.alpha(prop.par[j]), mean = mean.alpha[j], sd = sd.alpha[j], log = TRUE)}), sapply(5:13, function(j) {dnorm(transfo.beta(prop.par[j]), mean = mean.beta[j-4], sd = sd.beta[j-4], log=TRUE)}) )
		    }    
		    if(length(prop.par)==39){
		        vect.res = c( sapply(1:11, function(j) {dnorm(transfo.alpha(prop.par[j]), mean = mean.alpha[j], sd = sd.alpha[j], log = TRUE)}), sapply(12:39, function(j) {dnorm(transfo.beta(prop.par[j]), mean = mean.beta[j-11], sd = sd.beta[j-11], log=TRUE)}) )
		    }
	        return(ifelse(log, sum(vect.res), exp(sum(vect.res))))
	    }
	    stop("wrong type specification")
	}


	if(model=="Pairwise"){return( proposal.pb(type, cur.par, prop.par, MCpar, log) )}
	if(model=="Husler"){return( proposal.hr(type, cur.par, prop.par, MCpar, log) )}
	if(model=="Dirichlet"){return( proposal.di(type, cur.par, prop.par, MCpar, log) )}
	if(model=="Extremalt"){return( proposal.et(type, cur.par, prop.par, MCpar, log) )}
	if(model=="Asymmetric"){return( proposal.al(type, cur.par, prop.par, MCpar, log) )}				

}




posteriorMCMC <- function (
	Nsim, Nbin = 0, 
	Hpar, MCpar, 
	dat, par.start = NULL, 
	show.progress = floor(seq(1,Nsim, length.out = 20)), 
	seed = NULL, kind = "Mersenne-Twister", 
    save = FALSE, name.save = NULL, save.directory = "~", 
    name.dat = "", 
    model, c=NULL) 
{

	# Preliminary function
	
	lAccept.ratio <- function (cur.par, prop.par, llh.cur, lprior.cur, Hpar, MCpar, dat, model,c) {
 
	    p <- ncol(dat)
   
		new.ll <- dens(x = dat, model=model, par = prop.par, c=c, log = TRUE, vectorial = FALSE)
		        
	    new.lprior <- prior(model,type = "d", par = prop.par, Hpar = Hpar, log = TRUE, dimData = p)
	    proposal.oldToProp <- proposal(model,type = "d", cur.par = cur.par,  prop.par = prop.par, MCpar = MCpar, log = TRUE)
    	proposal.propToOld <- proposal(model,type = "d", cur.par = prop.par, 
		prop.par = cur.par, MCpar = MCpar, log = TRUE)
    	return(list(lrho = new.ll - llh.cur + new.lprior - lprior.cur + proposal.propToOld - proposal.oldToProp, llh = new.ll, 
        	lprior = new.lprior))
	}
	
	#

	argnames <- ls()
    arglist <- list()
    for (i in 1:length(argnames)) {
        arglist[[i]] <- get(argnames[i])
    }
    names(arglist) <- argnames
    if (!is.null(seed)) 
        set.seed(seed, kind = kind)
    start.time <- proc.time()
    p <- ncol(dat)
    if (is.null(par.start)) {
        while(!is.finite(condit)){
            par.start <- prior(model, type = "r", n = 1, Hpar = Hpar, dimData = p)
			condit1 <- dens(x = dat, model=model, par = par.start, c=c, log = TRUE,  vectorial = FALSE)
            condit <- condit1 + prior(model, type = "d", par = par.start, Hpar = Hpar, log = TRUE, dimData = p)
        }
    }
    cur.par <- par.start
    
	llh.cur <- dens(x = dat, model=model, par = cur.par, c=c, log = TRUE, vectorial = FALSE)
    
    lprior.cur <- prior(model, type = "d", par = cur.par, Hpar = Hpar, log = TRUE, dimData = p)
    nsim = 1
    n.accept = 0
    n.accept.kept = 0
    leng = length(cur.par)
    mean.res = rep(0, leng)
    emp.variance = rep(0, leng)
    emp.variance.unNorm = rep(0, leng)
    stored.vals <- matrix(0, nrow = Nsim - Nbin, ncol = leng)
    llh <- double(Nsim - Nbin)
    lprior <- double(Nsim - Nbin)
    stored.vals[1, ] = cur.par
    lprior[1] <- lprior.cur
    llh[1] <- llh.cur
    while (nsim <= Nsim) {
        if (any(nsim == show.progress)) {
            cat(paste("iter", nsim, ": n.accepted=", n.accept, "\n", sep = " "))
        }
        prop.par <- proposal(model,type = "r", cur.par = cur.par, prop.par = NULL, MCpar = MCpar)
        ratio.list <- lAccept.ratio(cur.par = cur.par, prop.par = prop.par, llh.cur = llh.cur, lprior.cur = lprior.cur, 
        	Hpar = Hpar, MCpar= MCpar, dat = dat, model = model)
        if ((is.finite(ratio.list$lrho)) && ((ratio.list$lrho > 0) || (log(runif(1)) <= ratio.list$lrho))) {
            n.accept <- n.accept + 1
            if (nsim > Nbin) 
                n.accept.kept = n.accept.kept + 1
            cur.par <- prop.par
            llh.cur <- ratio.list$llh
            lprior.cur <- ratio.list$lprior
        }
        if (nsim > Nbin) {
            n <- nsim - Nbin
            new.mean.res = mean.res + 1/n * (cur.par - mean.res)
            emp.variance.unNorm <- emp.variance.unNorm + (cur.par - new.mean.res) * (cur.par - mean.res)
            mean.res <- new.mean.res
            if (nsim == Nsim) {
                emp.variance <- emp.variance.unNorm/(n - 1)
            }
            stored.vals[n, ] <- cur.par
            llh[n] <- llh.cur
            lprior[n] <- lprior.cur
        }
        nsim <- nsim + 1
    }
    end.time <- proc.time()
    print(end.time - start.time)
    
	BIC <- -2*dens(x=dat, model=model, par=mean.res, c=c, log=TRUE, vectorial=FALSE)+length(mean.res)*(log(nrow(dat))+log(2*pi))
	   
    res <- list(stored.vals = stored.vals, llh = llh, lprior = lprior, 
        arguments = arglist, elapsed = end.time - start.time, 
        Nsim = Nsim, Nbin = Nbin, n.accept = n.accept, n.accept.kept = n.accept.kept, 
        emp.mean = mean.res, emp.sd = sqrt(emp.variance), BIC = BIC)
	class (res) <- "postsample"
	if (save && !is.null(name.save)) {
        loglist <- list(model = model, Nsim = Nsim, Nbin = Nbin, 
            dat = name.dat, Hpar = Hpar, MCpar = MCpar, seed = seed, 
            kind = kind)
        name.log <- paste(name.save, ".log", sep = "")
        assign(name.save, res)
        assign(name.log, loglist)
        save(list = name.save, file = paste(save.directory, "/", name.save, ".rda", sep = ""))
        save(list = name.log, file = paste(save.directory, "/", name.log, ".rda", sep = ""))
    }
    return(res)
}


