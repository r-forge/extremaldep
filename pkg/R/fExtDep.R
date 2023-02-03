fExtDep <- function(method="PPP", data, model, par.start = NULL, c = 0,
                       optim.method = "BFGS", trace = 0, sig = 3,
                       Nsim, Nbin = 0, Hpar, MCpar, seed = NULL
                       ){
  
  # Checking the method
  
  methods <- c("PPP", "BayesianPPP", "Composite")
  if(!any(method ==  methods)){ stop("estimation method wrongly specified")}
  
  # Checking the model
  
  if(method == "PPP" || method == "BayesianPPP"){
      models <- c("PB", "HR", "ET", "EST", "TD", "AL")
      if(!any(model ==  models)){ stop("model wrongly specified")}
  }else if(method == "Composite"){
      models <- c("HR", "ET", "EST")
      if(!any(model ==  models)){ stop("model wrongly specified")}
  }
  
  if(method == "PPP"){
  
    if(model == "AL" && c==0 && length(par.start)==4){ # 3-dimensional AL model with only parameters for the interior
      optimfun <- function(para){
        return(dExtDep(x=data, method="Parametric", model=model, par=c(rep(1,3), para[1], rep(0,6), para[2:4]), angular=TRUE, log=TRUE, c=c, vectorial=FALSE))
      }
    }else{
      optimfun <- function(para){
        return(dExtDep(x=data, method="Parametric", model=model, par=para, angular=TRUE, log=TRUE, c=c, vectorial=FALSE))
      }  
    }
    
    est.para <- optim(par.start, optimfun, method=optim.method, control=list(maxit=1e8, trace=trace, fnscale=-1), hessian=TRUE) 
    
    LogLik <- est.para$value # log likelihood
    param.est <- est.para$par # Parameters estimates
    J <- -est.para$hessian

    n <- nrow(data)
    #s <- 0
    score <- matrix(nrow=n,ncol=length(param.est))
    for(i in 1:n){
      stepJ <- function(para){
        if(model == "AL" && c==0 && length(par.start)==4){
          dExtDep(x=data[i,], method="Parametric", model=model, par=c(rep(1,3), para[1], rep(0,6), para[2:4]), angular=TRUE, log=TRUE, c=c, vectorial=FALSE)
        }else{
          dExtDep(x=data[i,], method="Parametric", model=model, par=para, angular=TRUE, log=TRUE, c=c, vectorial=FALSE)  
        }
      }
      #score[i,] = numDeriv::jacobian(stepJ,param.est)
      score[i,] =  fdjacobian(fun=stepJ, par=param.est, split=FALSE, epsVec=rep(0.005, length(param.est)))
      #s = s + numDeriv::hessian(stepJ,param.est)
    } 
    
    Ind <- apply(score, 1, function(x) any(x==0) || any(is.na(x))) # Indicate those wit 0s or NAs. Can happen when c !=0 and there are issues with the integral.
    K=var(score[!Ind,]); # variability matrix
    #J=-s; # sensitivity matrix
    sJ = solve(J);
    TIC = 2*matrix.trace(K %*% sJ * n)-2*LogLik; # TIC
    #SE = diag(matrix.sqrt(solve((J %*% solve(K) %*% J)/n ))) # Standard errors
    SE = diag(matrix.sqrt(sJ %*% K %*% sJ * n )) # Standard errors
    
    return(list(par=round(param.est,sig), LL=round(LogLik,sig), 
                  TIC=round(TIC,sig), SE=round(SE,sig) ))

  }
  if(method == "BayesianPPP"){
    fit <- posteriorMCMC(Nsim = Nsim, Nbin = Nbin, Hpar = Hpar, MCpar = MCpar, 
                  dat = data, par.start = par.start, seed = seed, 
                  model = model, c = c)
    return(fit)
  }
  if( method == "Composite"){
    
    biv.cllik <- function(par, data, model){ 
      if(model=="HR"){
        fun <- function(x, par){ pair.cllik.hr(z=x, par=par) }
      }
      if(model=="ET"){
        fun <- function(x, par){ pair.cllik.ext(z=x, par=par) }
      }
      if(model=="EST"){
        fun <- function(x, par){ pair.cllik.skewt(z=x, par=par) }
      }

      if(is.matrix(data)){
        return( sum( apply(data, 1, function(x){ fun(x,par)} ) ) )        
      }else if(is.vector(data)){
        return( fun(data,par) )
      }      

      
    }
    
    est <- optim(par.start, biv.cllik, data=data, model=model, method=optim.method, hessian=TRUE, control = list(maxit = 1e8, trace = trace, fnscale=-1)) 	
    
    n <- nrow(data)
    #s <- 0
    J <- -est$hessian
    
    score <- matrix(nrow=n, ncol=length(est$par))
    for(i in 1:n){
      stepJ <- function(para){
        biv.cllik(par=para, data=data[i,], model=model)
      }
      #score[i,] = numDeriv::jacobian(stepJ, est$par)
      score[i,] =  fdjacobian(fun=stepJ, par=est$par, split=FALSE, epsVec=rep(0.005, length(est$par)))
      #s = s + numDeriv::hessian(stepJ, est$par)	
      
    } 
    K=var(score); # variability matrix
    #J=-s; # sensitivity matrix
    sJ = solve(J);
    TIC = 2*matrix.trace(K %*% sJ * n)-2*est$value; # TIC
    SE = diag(matrix.sqrt(sJ %*% K %*% sJ * n )) # Standard errors

    return(list(par=round(est$par, sig), LL=round(est$value,sig),
                SE=round(SE,sig), TIC=round(TIC,sig) ))      
  }  

}

###############################################################################
###############################################################################
## Hidden functions 
###############################################################################
###############################################################################

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

###############################################################################
###############################################################################
## Hidden functions when method == ""BayesianPPP"
###############################################################################
###############################################################################


prior <- function(model, type = c("r", "d"), n, par, Hpar, log, dimData, c){
  
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
      #rho <- sqrt( invlogit( matrix(rnorm((lengthPar-1)*n, mean = Hpar$mean.rho, sd = Hpar$sd.rho),ncol=lengthPar-1)))
      rho <- tanh( matrix(rnorm((lengthPar-1)*n, mean = Hpar$mean.rho, sd = Hpar$sd.rho),ncol=lengthPar-1))
      nu <- exp(rnorm(n, mean = Hpar$mean.nu, sd = Hpar$sd.nu))
      res <- cbind(rho,nu)
      return(res)
    }
    if (type == "d") {
      #lpar <- c( logit(par[-length(par)]^2), log(par[length(par)]) )
      lpar <- c( atanh(par[-length(par)]), log(par[length(par)]) )
      ld1 <- dnorm(lpar[-length(lpar)], mean = Hpar$mean.rho, sd = Hpar$sd.rho, log = TRUE )
      ld2 <- dnorm(lpar[length(lpar)], mean = Hpar$mean.nu, sd = Hpar$sd.nu, log = TRUE )
      if (log) 
        return( sum(ld1) + ld2 )
      else return( exp( sum(ld1) + ld2 ) )
    }
    stop("wrong 'type' argument")
  }
  
  ###
  ### Skew-t
  ###
  
  prior.skt <- function(type = c("r", "d"), n, par, Hpar, log, dimData){
    if (type == "r") {
      p <- dimData
      # rho is of length choose(p,2), alpha is of length p and mu of length 1
      rho <- tanh( matrix(rnorm(choose(p,2)*n, mean = Hpar$mean.rho, sd = Hpar$sd.rho),ncol=choose(p,2) ))
      #rho <- sqrt( invlogit( matrix(rnorm(choose(p,2)*n, mean = Hpar$mean.rho, sd = Hpar$sd.rho), ncol=choose(p,2) )))
      alpha <- matrix(rnorm(p*n, mean = Hpar$mean.alpha, sd = Hpar$sd.alpha), ncol=p )
      nu <- exp(rnorm(n, mean = Hpar$mean.nu, sd = Hpar$sd.nu))
      res <- cbind(rho,alpha,nu)
      return(res)
    }
    if (type == "d") {
      lpar <- c( atanh(par[1:choose(dimData,2)]), par[(choose(dimData,2)+1):(length(par)-1)], log(par[length(par)]) )
      #lpar <- c( logit(par[1:choose(dimData,2)]^2), par[(choose(dimData,2)+1):(length(par)-1)], log(par[length(par)]) )
      ld1 <- dnorm(lpar[1:choose(dimData,2)], mean = Hpar$mean.rho, sd = Hpar$sd.rho, log = TRUE )
      ld2 <- dnorm(lpar[(choose(dimData,2)+1):(length(par)-1)], mean = Hpar$mean.alpha, sd = Hpar$sd.alpha, log = TRUE )
      ld3 <- dnorm(lpar[length(lpar)], mean = Hpar$mean.nu, sd = Hpar$sd.nu, log = TRUE )
      if (log)
        return( sum(ld1) + sum(ld2) + ld3 )
      else return( exp( sum(ld1) + sum(ld2) + ld3 ) )
    }
    stop("wrong 'type' argument")
  }
  
  ###
  ### Asymmetric Logistic
  ###
  
  
  prior.al <- function(type = c("r", "d"), n, par, Hpar, log, dimData, c){
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
        if(c==0){
          lengthPar <- 4
          alpha <- exp(matrix(rnorm(n, mean = Hpar$mean.alpha, sd = Hpar$sd.alpha),ncol=1))+1
          beta <- invlogit( matrix(rnorm(3*n, mean = Hpar$mean.beta, sd = Hpar$sd.beta),ncol=3))
        }else if(c>0){
          lengthPar <- 13
          alpha <- exp(matrix(rnorm(4*n, mean = Hpar$mean.alpha, sd = Hpar$sd.alpha),ncol=4))+1
          beta <- invlogit( matrix(rnorm(9*n, mean = Hpar$mean.beta, sd = Hpar$sd.beta),ncol=9))
        }
        res <- cbind(alpha,beta)
        return(res)
      }
      if(p==4){
        if(c==0){
          lengthPar <- 5
          alpha <- exp(matrix(rnorm(n, mean = Hpar$mean.alpha, sd = Hpar$sd.alpha),ncol=1))+1
          beta <- invlogit( matrix(rnorm(4*n, mean = Hpar$mean.beta, sd = Hpar$sd.beta),ncol=4))          
        }else if(c>0){
          lengthPar <- 39
          alpha <- exp(matrix(rnorm(11*n, mean = Hpar$mean.alpha, sd = Hpar$sd.alpha),ncol=11))+1
          beta <- invlogit( matrix(rnorm(28*n, mean = Hpar$mean.beta, sd = Hpar$sd.beta),ncol=28))
        }
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
        if(c==0){
          lpar <- c( log( par[1]-1), logit(par[2:4]))
          ld1 <- dnorm(lpar[1], mean = Hpar$mean.alpha, sd = Hpar$sd.alpha, log = TRUE)
          ld2 <- dnorm(lpar[2:4], mean = Hpar$mean.beta, sd = Hpar$sd.beta, log = TRUE)
        }else if(c>0){
          lpar <- c( log( par[1:4]-1), logit(par[5:13]))
          ld1 <- dnorm(lpar[1:4], mean = Hpar$mean.alpha, sd = Hpar$sd.alpha, log = TRUE)
          ld2 <- dnorm(lpar[5:13], mean = Hpar$mean.beta, sd = Hpar$sd.beta, log = TRUE)
        }
      }
      if(p==4){
        if(c==0){
          lpar <- c( log( par[1]-1), logit(par[2:5]) )
          ld1 <- dnorm(lpar[1], mean = Hpar$mean.alpha, sd = Hpar$sd.alpha, log = TRUE)
          ld2 <- dnorm(lpar[2:5], mean = Hpar$mean.beta, sd = Hpar$sd.beta, log = TRUE)
          
        }else if(c>0){
          lpar <- c( log( par[1:11]-1), logit(par[12:39]) )
          ld1 <- dnorm(lpar[1:11], mean = Hpar$mean.alpha, sd = Hpar$sd.alpha, log = TRUE)
          ld2 <- dnorm(lpar[12:39], mean = Hpar$mean.beta, sd = Hpar$sd.beta, log = TRUE)
        }
      }        
      if (log) 
        return( sum(ld1) + sum(ld2) )
      else return( exp( sum(ld1) + sum(ld2) ) )
    }
    stop("wrong 'type' argument")
  }
  
  if(model=="PB"){return( prior.pb(type, n, par, Hpar, log, dimData) )}
  if(model=="HR"){return( prior.hr(type, n, par, Hpar, log, dimData) )}
  if(model=="TD"){return( prior.di(type, n, par, Hpar, log, dimData) )}
  if(model=="ET"){return( prior.et(type, n, par, Hpar, log, dimData) )}
  if(model=="EST"){return( prior.skt(type, n, par, Hpar, log, dimData) )}
  if(model=="AL"){return( prior.al(type, n, par, Hpar, log, dimData,c) )}
  
}


proposal <- function (model, type = c("r", "d"), cur.par, prop.par, MCpar, log = TRUE,c) { 
  
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
    
    # transfo.rho <- function(x) { logit(x^2) }
    # invtransfo.rho <- function(x) { sign(x)*sqrt(invlogit(x)) }
    transfo.rho <- function(x){atanh(x)}
    invtransfo.rho <- function(x){ tanh(x)}
    transfo.df <- function(x){log(x)}
    invtransfo.df <- function(x){exp(x)}
    
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
  
  proposal.skt <- function (type = c("r", "d"), cur.par, prop.par, MCpar, log = TRUE) {
    
    Dim <- function(x){
      round(uniroot( function(y) choose(y,2) + y + 1 -length(x), interval=c(2,5)  )$root)
    }
    p <- Dim(cur.par)
    
    sd.rho <- rep(MCpar, choose(p,2))
    sd.alpha <- rep(MCpar, p)
    sd.df <-  MCpar
    
    transfo.rho <- function(x){atanh(x)}
    invtransfo.rho <- function(x){ tanh(x)}
    # transfo.rho <- function(x) { logit(x^2) }
    # invtransfo.rho <- function(x) { sign(x)*sqrt(invlogit(x)) }
    transfo.df <- function(x) { log(x) }
    invtransfo.df <- function(x) { exp(x) }
    
    mean.rho <- transfo.rho(cur.par[1:choose(p,2)])
    mean.alpha <- cur.par[(choose(p,2)+1):(length(cur.par)-1)]
    mean.df <- transfo.df(cur.par[length(cur.par)])
    if (type == "r") {
      res <- c(invtransfo.rho(rnorm(choose(p,2), mean = mean.rho, sd = sd.rho)),
               rnorm(p, mean = mean.alpha, sd = sd.alpha),
               invtransfo.df(rnorm(1, mean = mean.df, sd = sd.df)))
      return(res)
    }
    if (type == "d") {
      vect.res = c(dnorm(transfo.rho(prop.par[1:choose(p,2)]), mean = mean.rho, sd = sd.rho, log=TRUE),
                   dnorm(prop.par[(choose(p,2)+1):(length(prop.par)-1)], mean = mean.alpha, sd = sd.alpha, log=TRUE),
                   dnorm(transfo.df(prop.par[length(prop.par)]), mean = mean.df, sd = sd.df, log = TRUE))
      return(ifelse(log, sum(vect.res), exp(sum(vect.res))))
    }
    stop("wrong type specification")
  }
  
  proposal.al <- function (type = c("r", "d"), cur.par, prop.par, MCpar, log = TRUE, c) {
    
    if(length(cur.par)==3){ #2-d
      sd.alpha <- rep(MCpar, 1)
      sd.beta <-  rep(MCpar, 2)
    }
    if(c==0){
      if(length(cur.par)==4){ #3-d
        sd.alpha <- rep(MCpar, 1)
        sd.beta <-  rep(MCpar, 3)
      }
      if(length(cur.par)==5){ #4-d
        sd.alpha <- rep(MCpar, 1)
        sd.beta <-  rep(MCpar, 4)
      } 
    }else if(c>0){
      if(length(cur.par)==13){ #3-d
        sd.alpha <- rep(MCpar, 4)
        sd.beta <-  rep(MCpar, 9)
      }
      if(length(cur.par)==39){ #4-d
        sd.alpha <- rep(MCpar, 11)
        sd.beta <-  rep(MCpar, 28)
      }    
    }
    
    
    transfo.alpha <- function(x) { log(x-1) }
    invtransfo.alpha <- function(x) { exp(x)+1 }    
    transfo.beta <- function(x) { logit(x) }
    invtransfo.beta <- function(x) { invlogit(x) }
    
    if(length(cur.par)==3){
      mean.alpha <- transfo.alpha(cur.par[1])
      mean.beta <- transfo.beta(cur.par[2:3])
    }
    if(c==0){
      if(length(cur.par)==4){
        mean.alpha <- transfo.alpha(cur.par[1])
        mean.beta <- transfo.beta(cur.par[2:4])
      }
      if(length(cur.par)==5){
        mean.alpha <- transfo.alpha(cur.par[1])
        mean.beta <- transfo.beta(cur.par[2:5])
      }	
    }else if(c>0){
      if(length(cur.par)==13){
        mean.alpha <- transfo.alpha(cur.par[1:4])
        mean.beta <- transfo.beta(cur.par[5:13])
      }
      if(length(cur.par)==39){
        mean.alpha <- transfo.alpha(cur.par[1:11])
        mean.beta <- transfo.beta(cur.par[12:39])
      }		
    }
    if (type == "r") {
      if(length(cur.par)==3){
        res <- c(invtransfo.alpha(rnorm(1, mean = mean.alpha, sd = sd.alpha)), invtransfo.beta(rnorm(2, mean = mean.beta, sd = sd.beta)))
      }	
      if(c==0){
        if(length(cur.par)==4){
          res <- c(invtransfo.alpha(rnorm(1, mean = mean.alpha, sd = sd.alpha)), invtransfo.beta(rnorm(3, mean = mean.beta, sd = sd.beta)))
        }
        if(length(cur.par)==5){
          res <- c(invtransfo.alpha(rnorm(1, mean = mean.alpha, sd = sd.alpha)), invtransfo.beta(rnorm(4, mean = mean.beta, sd = sd.beta)))
        }
      }else if(c>0){
        if(length(cur.par)==13){
          res <- c(invtransfo.alpha(rnorm(4, mean = mean.alpha, sd = sd.alpha)), invtransfo.beta(rnorm(9, mean = mean.beta, sd = sd.beta)))
        }
        if(length(cur.par)==39){
          res <- c(invtransfo.alpha(rnorm(11, mean = mean.alpha, sd = sd.alpha)), invtransfo.beta(rnorm(28, mean = mean.beta, sd = sd.beta)))
        }
      }
      
      return(res)        
    }
    if (type == "d") {
      if(length(prop.par)==3){
        vect.res = c( dnorm(transfo.alpha(prop.par[1]), mean = mean.alpha, sd = sd.alpha, log = TRUE), sapply(2:3, function(j) {dnorm(transfo.beta(prop.par[j]), mean = mean.beta[j-1], sd = sd.beta[j-1], log=TRUE)}) )
      }
      if(c==0){
        if(length(prop.par)==4){
          vect.res = c( dnorm(transfo.alpha(prop.par[1]), mean = mean.alpha[1], sd = sd.alpha[1], log = TRUE), sapply(2:4, function(j) {dnorm(transfo.beta(prop.par[j]), mean = mean.beta[j-1], sd = sd.beta[j-1], log=TRUE)}) )
        }    
        if(length(prop.par)==5){
          vect.res = c( dnorm(transfo.alpha(prop.par[1]), mean = mean.alpha[1], sd = sd.alpha[1], log = TRUE), sapply(2:5, function(j) {dnorm(transfo.beta(prop.par[j]), mean = mean.beta[j-1], sd = sd.beta[j-1], log=TRUE)}) )
        }
      }else if(c>0){
        if(length(prop.par)==13){
          vect.res = c( sapply(1:4, function(j) {dnorm(transfo.alpha(prop.par[j]), mean = mean.alpha[j], sd = sd.alpha[j], log = TRUE)}), sapply(5:13, function(j) {dnorm(transfo.beta(prop.par[j]), mean = mean.beta[j-4], sd = sd.beta[j-4], log=TRUE)}) )
        }    
        if(length(prop.par)==39){
          vect.res = c( sapply(1:11, function(j) {dnorm(transfo.alpha(prop.par[j]), mean = mean.alpha[j], sd = sd.alpha[j], log = TRUE)}), sapply(12:39, function(j) {dnorm(transfo.beta(prop.par[j]), mean = mean.beta[j-11], sd = sd.beta[j-11], log=TRUE)}) )
        }
      }
      return(ifelse(log, sum(vect.res), exp(sum(vect.res))))
    }
    stop("wrong type specification")
  }
  
  
  if(model=="PB"){return( proposal.pb(type, cur.par, prop.par, MCpar, log) )}
  if(model=="HR"){return( proposal.hr(type, cur.par, prop.par, MCpar, log) )}
  if(model=="TD"){return( proposal.di(type, cur.par, prop.par, MCpar, log) )}
  if(model=="ET"){return( proposal.et(type, cur.par, prop.par, MCpar, log) )}
  if(model=="EST"){return( proposal.skt(type, cur.par, prop.par, MCpar, log) )}
  if(model=="AL"){return( proposal.al(type, cur.par, prop.par, MCpar, log, c))}				
  
}




posteriorMCMC <- function (
    Nsim, Nbin = 0, 
    Hpar, MCpar, 
    dat, par.start = NULL,
    seed = NULL,
    model, c=0)
{
  
  # Preliminary function
  
  lAccept.ratio <- function (cur.par, prop.par, llh.cur, lprior.cur, Hpar, MCpar, dat, model,c) {
    
    p <- ncol(dat)
    
    ll.fun <- function(para, model, c){
      if(model == "AL" && c==0 && length(para)==4){ # 3-dimensional AL model with only parameters for the interior
        return(dExtDep(x=dat, method="Parametric", model=model, par=c(rep(1,3), para[1], rep(0,6), para[2:4]), angular=TRUE, log=TRUE, c=c, vectorial=FALSE))
      }else if(model == "AL" && c==0 && length(par.start)==5){ # 4-dimensional AL model with only parameters for the interior
        return(dExtDep(x=dat, method="Parametric", model=model, par=c(rep(1,10), para[1], rep(0,24), para[2:4]), angular=TRUE, log=TRUE, c=c, vectorial=FALSE))
      }else{
        return(dExtDep(x=dat, method="Parametric", model=model, par=para, angular=TRUE, log=TRUE, c=c, vectorial=FALSE))
      }
    }  
    
    #new.ll <- dExtDep(x = dat, method="Parametric", model=model, par=prop.par, angular=TRUE, log=TRUE, c=c, vectorial=FALSE)
    new.ll <- ll.fun(para=prop.par, model=model, c=c)
    
    new.lprior <- prior(model,type = "d", par = prop.par, Hpar = Hpar, log = TRUE, dimData = p, c = c)
    proposal.oldToProp <- proposal(model,type = "d", cur.par = cur.par,  prop.par = prop.par, MCpar = MCpar, log = TRUE, c = c)
    proposal.propToOld <- proposal(model,type = "d", cur.par = prop.par, prop.par = cur.par, MCpar = MCpar, log = TRUE, c = c)
    return(list(lrho = new.ll - llh.cur + new.lprior - lprior.cur + proposal.propToOld - proposal.oldToProp, llh = new.ll, 
                lprior = new.lprior))
  }
  
  ###
  
  ll.fun <- function(para, model, c){
    if(model == "AL" && c==0 && length(para)==4){ # 3-dimensional AL model with only parameters for the interior
        return(dExtDep(x=dat, method="Parametric", model=model, par=c(rep(1,3), para[1], rep(0,6), para[2:4]), angular=TRUE, log=TRUE, c=c, vectorial=FALSE))
    }else if(model == "AL" && c==0 && length(par.start)==5){ # 4-dimensional AL model with only parameters for the interior
        return(dExtDep(x=dat, method="Parametric", model=model, par=c(rep(1,10), para[1], rep(0,24), para[2:4]), angular=TRUE, log=TRUE, c=c, vectorial=FALSE))
    }else{
        return(dExtDep(x=dat, method="Parametric", model=model, par=para, angular=TRUE, log=TRUE, c=c, vectorial=FALSE))
    }
  }  
  
  ###
  
  argnames <- ls()
  arglist <- list()
  for (i in 1:length(argnames)) {
    arglist[[i]] <- get(argnames[i])
  }
  names(arglist) <- argnames
  if (!is.null(seed)) 
    set.seed(seed)
  start.time <- proc.time()
  p <- ncol(dat)
  if (is.null(par.start)) {
    repeat{
      par.start <- prior(model, type = "r", n = 1, Hpar = Hpar, dimData = p, c = c)
      #condit <- dExtDep(x=dat, method="Parametric", model=model, par=par.start, angular=TRUE, log=TRUE, c=c, vectorial=FALSE)
      condit <- ll.fun(para=par.start, model=model, c=c)
      if(is.finite(condit)) break
    }
  }
  cur.par <- par.start
  
  #llh.cur <- dExtDep(x=dat, method="Parametric", model=model, par=cur.par, angular=TRUE, log=TRUE, c=c, vectorial=FALSE)
  llh.cur <- ll.fun(para=cur.par, model=model, c=c)
  lprior.cur <- prior(model, type = "d", par = cur.par, Hpar = Hpar, log = TRUE, dimData = p, c = c)  
  
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
  
  pb = txtProgressBar(min = 0, max = Nsim, initial = 0, style=3) #creates a progress bar
  print.nsim <- seq(0, Nsim, by=100)
  
  while (nsim <= Nsim) {
    
    prop.par <- proposal(model,type = "r", cur.par = cur.par, prop.par = NULL, MCpar = MCpar, c = c)
    ratio.list <- lAccept.ratio(cur.par = cur.par, prop.par = prop.par, llh.cur = llh.cur, lprior.cur = lprior.cur, 
                                Hpar = Hpar, MCpar= MCpar, dat = dat, model = model, c=c)
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
    
    if(nsim %in% print.nsim) setTxtProgressBar(pb, nsim)
  }
  end.time <- proc.time()
  
  close(pb)
  
  #BIC <- -2*dExtDep(x=dat, method="Parametric", model=model, par=mean.res, angular=TRUE, c=c, log=TRUE, vectorial=FALSE)+length(mean.res)*(log(nrow(dat))+log(2*pi))
  BIC <- -2*ll.fun(para=mean.res, model=model, c=c)+length(mean.res)*(log(nrow(dat))+log(2*pi))
  
  res <- list(stored.vals = stored.vals, llh = llh, lprior = lprior, 
              arguments = arglist, elapsed = end.time - start.time, 
              Nsim = Nsim, Nbin = Nbin, n.accept = n.accept, n.accept.kept = n.accept.kept, 
              emp.mean = mean.res, emp.sd = sqrt(emp.variance), BIC = BIC)
  class (res) <- "postsample"
  
  return(res)
}

###############################################################################
###############################################################################
## Hidden functions when method == "Composite"
###############################################################################
###############################################################################

###############################################################################
## Routines for pairwise likelihood for Husler Reiss model 
###############################################################################

llHR <- function(x, lambda=1){
  if(is.matrix(x)){
    n <- nrow(x)
  }else{ n <- 1}
  par <- c(df,scale)
  .C("llHRmax", as.double(x), as.double(lambda), as.integer(n), out=double(1), NAOK=TRUE)$out
  
}

pair.cllik.hr <- function(z,par){
  
  d <- length(z)
  pairs <- t(combn(d,2))
  
  pllik <- sum(sapply(1:nrow(pairs), function(p) llHR(z[pairs[p,]], lambda = par[p]) ))
    
  return(pllik)
}

###############################################################################
## Routines for pairwise likelihood for Extremal-t model 
###############################################################################

llET <- function(x, scale=1, df=1){
  if(is.matrix(x)){
    n <- nrow(x)
  }else{ n <- 1}
  par <- c(df,scale)
  .C("llETmax", as.double(x), as.double(par), as.integer(n), out=double(1), NAOK=TRUE)$out
}

pair.cllik.ext <- function(z,par){
  
  d <- length(z)
  pairs <- t(combn(d,2))
  
  pllik <- sum(sapply(1:nrow(pairs), function(p) llET(z[pairs[p,]], scale = par[p], df= tail(par,1) ) ))
  
  return(pllik)
}

###############################################################################
## Routines for pairwise likelihood for Extremal Skew-t model 
###############################################################################

llextst <- function(x, scale=1, shape=rep(0,length(x)), df=1){
  if(is.matrix(x)){
    n <- nrow(x)
  }else{ n <- 1}
  
  .C("llextst", as.double(x), as.integer(n), as.double(scale), as.double(df),
     as.double(shape), out=double(1), NAOK=TRUE)$out
}

pair.cllik.skewt <- function(z,par){
  
  d <- length(z)
  pairs <- t(combn(d,2))
  
  pllik <- sum(sapply(1:nrow(pairs), function(p) llextst(z[pairs[p,]], scale = par[p], shape = par[nrow(pairs)+pairs[p,]], df= tail(par,1) ) ))
  
  return(pllik)
}

