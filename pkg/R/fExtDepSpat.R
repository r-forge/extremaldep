fExtDepSpat <- function(model, z, sites, hit, jw, thresh, DoF, range, smooth, alpha, par0, 
                        acov1, acov2,
                        parallel, ncores, args1, args2, seed=123, 
                        method = "BFGS", sandwich=TRUE,
                        control = list(trace=1, maxit=50, REPORT=1, reltol=0.0001)){
  
  # models <- c("GG", "BR", "ET", "EST")
  # Currently limited to "ET" and "EST" !!!
  models <- c("ET", "EST")
  if(!any(model ==  models)){ stop("model wrongly specified")}
  
  Ns <- nrow(sites)
  if(missing(jw)){jw <- Ns} # Full likelihood
  if(jw > Ns){stop("'jw' cannot be greater than the number of 'sites'.")}
  if(jw < Ns){
    
    if(missing(thresh)){
      thresh <- max(dist(sites))+1
      message("'thresh' not provided, largest distance considered.")
    }
    
    cmat <- getjwise(sites = sites, jw=jw, thresh=thresh)
  }
  
  if(model == "ET"){
    
    if(missing(DoF) && missing(smooth) && missing(range)){
      if(length(par0) != 3){stop("'par0' should be of length 3")}
      par0 <- c(log(par0[1]), log(par0[2]), logit(par0[3]/2))
      control$parscale <- c(1,1,1)
    }
    
    if(missing(DoF) && missing(smooth) && !missing(range)){
      if(length(par0) != 2){stop("'par0' should be of length 2")}
      par0 <- c(log(par0[1]), logit(par0[2]/2))
      control$parscale <- c(1,1)
    }
    
    if(missing(DoF) && !missing(smooth) && missing(range)){
      if(length(par0) != 2){stop("'par0' should be of length 2")}
      par0 <- c(log(par0[1]), log(par0[2]))
      control$parscale <- c(1,1)
    }
    
    if(!missing(DoF) && missing(smooth) && missing(range)){
      if(length(par0) != 2){stop("'par0' should be of length 2")}
      par0 <- c(log(par0[1]), logit(par0[2]/2))
      control$parscale <- c(1,1)
    }
    
    if(!missing(DoF) && !missing(smooth) && missing(range)){
      if(length(par0) != 1){stop("'par0' should be of length 1")}
      par0 <- log(par0)
      control$parscale <- 1
    }
    
    if(!missing(DoF) && missing(smooth) && !missing(range)){
      if(length(par0) != 1){stop("'par0' should be of length 1")}
      par0 <- logit(par0/2)
      control$parscale <- 1
    }
    
    if(missing(DoF) && !missing(smooth) && !missing(range)){
      if(length(par0) != 1){stop("'par0' should be of length 1")}
      par0 <- log(par0)
      control$parscale <- 1
    }
    
  } # END if model == "ET"
  
  if(model == "EST"){
    
    if(missing(DoF) && missing(range) && missing(smooth) && missing(alpha) ){
      if(length(par0) != 6){stop("'par0' should be of length 6")}
      par0 <- c(log(par0[1]), log(par0[2]), logit(par0[3]/2), par0[4:6])
      control$parscale <- rep(1,6)
    }else if(!missing(DoF) && missing(range) && missing(smooth) && missing(alpha)){
      if(length(par0) != 5){stop("'par0' should be of length 5")}
      par0 <- c(log(par0[1]), logit(par0[2]/2), par0[3:5])
      control$parscale <- rep(1,5)
    }else if(missing(DoF) && !missing(range) && missing(smooth) && missing(alpha)){
      if(length(par0) != 5){stop("'par0' should be of length 5")}
      par0 <- c(log(par0[1]), logit(par0[2]/2), par0[3:5])
      control$parscale <- rep(1,5)
    }else if(missing(DoF) && missing(range) && !missing(smooth) && missing(alpha)){
      if(length(par0) != 5){stop("'par0' should be of length 5")}
      par0 <- c(log(par0[1]), log(par0[2]), par0[3:5])
      control$parscale <- rep(1,5)
    }else if(missing(DoF) && missing(range) && missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(length(par0) != 5){stop("'par0' should be of length 5")}
          par0 <- c(log(par0[1]), log(par0[2]), logit(par0[3]/2), par0[4:5])
          control$parscale <- rep(1,5)
        }else if(length(ind)==2){
          if(length(par0) != 4){stop("'par0' should be of length 4")}
          par0 <- c(log(par0[1]), log(par0[2]), logit(par0[3]/2), par0[4])
          control$parscale <- rep(1,4)
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        if(length(par0) != 3){stop("'par0' should be of length 3")}
        par0 <- c(log(par0[1]), log(par0[2]), logit(par0[3]/2))
        control$parscale <- rep(1,3)
      }
      
    }else if(!missing(DoF) && !missing(range) && missing(smooth) && missing(alpha)){
      if(length(par0) != 4){stop("'par0' should be of length 4")}
      par0 <- c(logit(par0[1]/2), par0[2:4])
      control$parscale <- rep(1,4)
    }else if(!missing(DoF) && missing(range) && !missing(smooth) && missing(alpha)){
      if(length(par0) != 4){stop("'par0' should be of length 4")}
      par0 <- c(log(par0[1]), par0[2:4])
      control$parscale <- rep(1,4)
    }else if(!missing(DoF) && missing(range) && missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(length(par0) != 4){stop("'par0' should be of length 4")}
          par0 <- c(log(par0[1]), logit(par0[2]/2), par0[3:4])
          control$parscale <- rep(1,4) 
        }else if(length(ind)==2){
          if(length(par0) != 3){stop("'par0' should be of length 3")}
          par0 <- c(log(par0[1]), logit(par0[2]/2), par0[3])
          control$parscale <- rep(1,3)
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        if(length(par0) != 2){stop("'par0' should be of length 2")}
        par0 <- c(log(par0[1]),logit(par0[2]/2))
        control$parscale <- rep(1,2)
      }
      
    }else if(missing(DoF) && !missing(range) && !missing(smooth) && missing(alpha)){
      if(length(par0) != 4){stop("'par0' should be of length 4")}
      par0 <- c(log(par0[1]), par0[3:4])
      control$parscale <- rep(1,4)
    }else if(missing(DoF) && !missing(range) && missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(length(par0) != 4){stop("'par0' should be of length 4")}
          par0 <- c(log(par0[1]),logit(par0[2]/2), par0[3:4])
          control$parscale <- rep(1,4)
        }else if(length(ind)==2){
          if(length(par0) != 3){stop("'par0' should be of length 3")}
          par0 <- c(log(par0[1]), logit(par0[2]/2), par0[3])
          control$parscale <- rep(1,3)
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        if(length(par0) != 2){stop("'par0' should be of length 2")}
        par0 <- c(log(par0[1]), logit(par0[2]/2))
        control$parscale <- rep(1,2)
      }
      
    }else if(missing(DoF) && missing(range) && !missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(length(par0) != 4){stop("'par0' should be of length 4")}
          par0 <- c(log(par0[1]), log(par0[2]), par0[3:4])
          control$parscale <- rep(1,4)
        }else if(length(ind)==2){
          if(length(par0) != 3){stop("'par0' should be of length 3")}
          par0 <- c(log(par0[1]), log(par0[2]), par0[3])
          control$parscale <- rep(1,3)
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        if(length(par0) != 2){stop("'par0' should be of length 2")}
        par0 <- c(log(par0[1]), log(par0[2]))
        control$parscale <- rep(1,2)
      }
      
    }else if(!missing(DoF) && !missing(range) && !missing(smooth) && missing(alpha)){
      if(length(par0) != 3){stop("'par0' should be of length 3")}
      par0 <- c(par0[1:3])
      control$parscale <- rep(1,3)
    }else if(!missing(DoF) && !missing(range) && missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(length(par0) != 3){stop("'par0' should be of length 3")}
          par0 <- c(logit(par0[1]/2), par0[2:3])
          control$parscale <- rep(1,3)
        }else if(length(ind)==2){
          if(length(par0) != 2){stop("'par0' should be of length 2")}
          par0 <- c(logit(par0[1]/2), par0[2])
          control$parscale <- rep(1,2)
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        if(length(par0) != 1){stop("'par0' should be of length 1")}
        par0 <- logit(par0/2)
        control$parscale <- rep(1,1)
      }
      
    }else if(!missing(DoF) && missing(range) && !missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(length(par0) != 3){stop("'par0' should be of length 3")}
          par0 <- c(log(par0[1]), par0[2:3])
          control$parscale <- rep(1,3)
        }else if(length(ind)==2){
          if(length(par0) != 2){stop("'par0' should be of length 2")}
          par0 <- c(log(par0[1]), par0[2])
          control$parscale <- rep(1,2)
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        if(length(par0) != 1){stop("'par0' should be of length 1")}
        par0 <- log(par0)
        control$parscale <- rep(1,1)
      }
      
    }else if(missing(DoF) && !missing(range) && !missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(length(par0) != 3){stop("'par0' should be of length 3")}
          par0 <- c(log(par0[1]), par0[2:3])
          control$parscale <- rep(1,3) 
        }else if(length(ind)==2){
          if(length(par0) != 2){stop("'par0' should be of length 2")}
          par0 <- c(log(par0[1]), par0[2])
          control$parscale <- rep(1,2)
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        if(length(par0) != 1){stop("'par0' should be of length 1")}
        par0 <- log(par0)
        control$parscale <- rep(1,1)
      }
      
    }else if(!missing(DoF) && !missing(range) && !missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(length(par0) != 2){stop("'par0' should be of length 2")}
          control$parscale <- rep(1,2) 
        }else if(length(ind)==2){
          if(length(par0) != 1){stop("'par0' should be of length 1")}
          control$parscale <- 1
        }
      }else{ 
        stop("All parameters can't be fixed")
      }
      
    }  
    
  } # END if model == "EST"
  
  
    
  ### Parallelisation
  if(parallel){
    if(missing(ncores)){
      message("'ncores' must be specified, setting it to 4 by default")
      ncores <- 4
    }
    
    cl <- makeCluster(ncores, type="PSOCK")
    registerDoParallel(cl)
  }
  
  if(model == "ET"){
    
    if(jw == Ns){ # Full likelihood
      
      if(missing(DoF) && missing(smooth) && missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>500){return(Inf)}
          if(exp(param[2])>50){return(Inf)}
          if(round(inv.logit(param[3]),6)==0 || round(inv.logit(param[3]),6)==1){return(Inf)}
          nllh_et(par=param, z = z, sites = sites, hit = hit, 
                  split = FALSE, parallel = parallel, 
                  pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(missing(DoF) && missing(smooth) && !missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>500){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          nllh_etr(par=param, range = range, z = z, sites = sites, hit = hit, 
                   split = FALSE, parallel = parallel, 
                   pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(missing(DoF) && !missing(smooth) && missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>500){return(Inf)}
          if(exp(param[2])>50){return(Inf)}
          nllh_ets(par=param, z = z, sites = sites, hit = hit, 
                   smooth = smooth, split = FALSE, parallel = parallel, 
                   pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(!missing(DoF) && missing(smooth) && missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[2])>50){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          nllh_etf(par=param, lDoF = log(DoF), z = z, sites = sites, hit = hit, 
                   split = FALSE, parallel = parallel, 
                   pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(!missing(DoF) && !missing(smooth) && missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param)>500){return(Inf)}
          nllh_etsf(par=param, lDoF=log(DoF), z = z, sites = sites, hit = hit, 
                    smooth = smooth, split = FALSE, parallel = parallel, 
                    pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(!missing(DoF) && missing(smooth) && !missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(round(inv.logit(param),6)==0 || round(inv.logit(param),6)==1){return(Inf)}
          nllh_etfr(par=param, lDoF=log(DoF), z = z, sites = sites, hit = hit, 
                    range = range, split = FALSE, parallel = parallel, 
                    pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(missing(DoF) && !missing(smooth) && !missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param)>500){return(Inf)}
          nllh_etsr(par=param, smooth = smooth, z = z, sites = sites, hit = hit, 
                    range = range, split = FALSE, parallel = parallel, 
                    pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }
      
    }else{ # j-wise composite likelihood
      
      if(missing(DoF) && missing(smooth) && missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>500){return(Inf)}
          if(exp(param[2])>50){return(Inf)}
          if(round(inv.logit(param[3]),6)==0 || round(inv.logit(param[3]),6)==1){return(Inf)}
          
          nllh_etjw(par=param, z = z, sites = sites, hit = hit,
                    jw = jw, cmat = cmat, 
                    split = FALSE, parallel = parallel, 
                    pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(missing(DoF) && missing(smooth) && !missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>500){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          nllh_etrjw(par=param, range = range, z = z, sites = sites, hit = hit, 
                     jw = jw, cmat = cmat,
                     split = FALSE, parallel = parallel, 
                     pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(missing(DoF) && !missing(smooth) && missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>500){return(Inf)}
          if(exp(param[2])>50){return(Inf)}
          nllh_etsjw(par=param, z = z, sites = sites, hit = hit, 
                     jw = jw, cmat = cmat,
                     smooth = smooth, split = FALSE, parallel = parallel, 
                     pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(!missing(DoF) && missing(smooth) && missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          nllh_etfjw(par=param, lDoF = log(DoF), z = z, sites = sites, hit = hit, 
                     jw = jw, cmat = cmat,
                     split = FALSE, parallel = parallel, 
                     pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(!missing(DoF) && !missing(smooth) && missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param)>500){return(Inf)}
          nllh_etsfjw(par=param, lDoF=log(DoF), z = z, sites = sites, hit = hit, 
                      jw = jw, cmat = cmat,
                      smooth = smooth, split = FALSE, parallel = parallel, 
                      pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(!missing(DoF) && missing(smooth) && !missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(round(inv.logit(param),6)==0 || round(inv.logit(param),6)==1){return(Inf)}
          nllh_etfrjw(par=param, lDoF=log(DoF), z = z, sites = sites, hit = hit, 
                      jw = jw, cmat = cmat,
                      range = range, split = FALSE, parallel = parallel, 
                      pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }else if(missing(DoF) && !missing(smooth) && !missing(range)){
        
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param)>50){return(Inf)}
          nllh_etsrjw(par=param, z = z, sites = sites, hit = hit, 
                      jw = jw, cmat = cmat,
                      smooth = smooth, range = range, split = FALSE, parallel = parallel, 
                      pfun = mypmvt, args1 = args1, args2 = args2, seed = seed) 
        }
        
      }
      
    } # END else jw != Ns  
    
  } # END if model == "ET"
  
  if(model == "EST"){
    
    if(jw == Ns){ # Full likelihood
      
      if(missing(DoF) && missing(range) && missing(smooth) && missing(alpha) ){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(exp(param[2])>500){return(Inf)}
          if(round(inv.logit(param[3]),6)==0 || round(inv.logit(param[3]),6)==1){return(Inf)}
          nllh_est(par=param, z = z, sites = sites, hit = hit, 
                   acov1=acov1, acov2=acov2,
                   split = FALSE, parallel = parallel, 
                   pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(!missing(DoF) && missing(range) && missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>500){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          nllh_estf(par=param, z = z, lDoF = log(DoF), sites = sites, hit = hit, 
                    acov1=acov1, acov2=acov2,
                    split = FALSE, parallel = parallel, 
                    pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(missing(DoF) && !missing(range) && missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          nllh_estr(par=param, z = z, range = range, sites = sites, hit = hit, 
                    acov1=acov1, acov2=acov2,
                    split = FALSE, parallel = parallel, 
                    pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(missing(DoF) && missing(range) && !missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(exp(param[2])>500){return(Inf)}
          nllh_ests(par=param, z = z, smooth = smooth, sites = sites, hit = hit, 
                    acov1=acov1, acov2=acov2,
                    split = FALSE, parallel = parallel, 
                    pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(missing(DoF) && missing(range) && missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(exp(param[2])>500){return(Inf)}
          if(round(inv.logit(param[3]),6)==0 || round(inv.logit(param[3]),6)==1){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_esta0(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                           acov1=acov1, acov2=acov2,
                           split = FALSE, parallel = parallel, 
                           pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_esta1(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                           acov1=acov1, acov2=acov2,
                           split = FALSE, parallel = parallel, 
                           pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_esta2(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                           acov1=acov1, acov2=acov2,
                           split = FALSE, parallel = parallel, 
                           pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_esta01(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                            acov1=acov1, acov2=acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_esta02(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                            acov1=acov1, acov2=acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_esta12(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                            acov1=acov1, acov2=acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_esta012(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                         acov1=acov1, acov2=acov2,
                         split = FALSE, parallel = parallel, 
                         pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
        }
      }else if(!missing(DoF) && !missing(range) && missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(round(inv.logit(param[1]),6)==0 || round(inv.logit(param[1]),6)==1){return(Inf)}
          nllh_estfr(par=param, z = z, lDoF = log(DoF), range = range, sites = sites, hit = hit, 
                    acov1=acov1, acov2=acov2,
                    split = FALSE, parallel = parallel, 
                    pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(!missing(DoF) && missing(range) && !missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>500){return(Inf)}
          nllh_estfs(par=param, z = z, lDoF = log(DoF), smooth = smooth, sites = sites, hit = hit, 
                    acov1 = acov1, acov2 = acov2,
                    split = FALSE, parallel = parallel, 
                    pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(!missing(DoF) && missing(range) && missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estfa0(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estfa1(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estfa2(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estfa01(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estfa02(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estfa12(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estfa012(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                         acov1 = acov1, acov2 = acov2,
                         split = FALSE, parallel = parallel, 
                         pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
        }
      }else if(missing(DoF) && !missing(range) && !missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          nllh_estrs(par=param, z = z, range = range, smooth = smooth, sites = sites, hit = hit, 
                     acov1 = acov1, acov2 = acov2,
                     split = FALSE, parallel = parallel, 
                     pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(missing(DoF) && !missing(range) && missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estra0(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estra1(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estra2(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estra01(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estra02(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estra12(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estra012(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                          acov1 = acov1, acov2 = acov2,
                          split = FALSE, parallel = parallel, 
                          pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
          
        }
      }else if(missing(DoF) && missing(range) && !missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(exp(param[2])>500){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estsa0(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estsa1(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estsa2(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estsa01(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estsa02(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estsa12(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estsa012(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                          acov1 = acov1, acov2 = acov2,
                          split = FALSE, parallel = parallel, 
                          pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
          
        }
      }else if(!missing(DoF) && !missing(range) && !missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          nllh_estfrs(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, sites = sites, hit = hit, 
                      acov1 = acov1, acov2 = acov2,
                      split = FALSE, parallel = parallel, 
                      pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(!missing(DoF) && !missing(range) && missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(round(inv.logit(param),6)==0 || round(inv.logit(param),6)==1){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estfra0(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estfra1(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estfra2(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estfra01(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estfra02(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estfra12(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estfra012(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                          acov1 = acov1, acov2 = acov2,
                          split = FALSE, parallel = parallel, 
                          pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
        }
      }else if(!missing(DoF) && missing(range) && !missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param)>500){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estfsa0(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estfsa1(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estfsa2(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estfsa01(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estfsa02(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estfsa12(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estfsa012(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                           acov1 = acov1, acov2 = acov2,
                           split = FALSE, parallel = parallel, 
                           pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
        }
      }else if(missing(DoF) && !missing(range) && !missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param)>50){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estrsa0(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estrsa1(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estrsa2(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estrsa01(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estrsa02(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estrsa12(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estrsa012(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                           acov1 = acov1, acov2 = acov2,
                           split = FALSE, parallel = parallel, 
                           pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
        }
      }else if(!missing(DoF) && !missing(range) && !missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estfrsa0(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estfrsa1(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estfrsa2(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estfrsa01(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estfrsa02(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estfrsa12(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            } # END if length(ind)==2
          } # END if any(is.na(alpha))
          
        } # END if nllh_tmp
      } 
      
    }else{ # Composite likelihood
      
      if(missing(DoF) && missing(range) && missing(smooth) && missing(alpha) ){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(exp(param[2])>500){return(Inf)}
          if(round(inv.logit(param[3]),6)==0 || round(inv.logit(param[3]),6)==1){return(Inf)}
          nllh_estjw(par=param, z = z, sites = sites, hit = hit, 
                     acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                     split = FALSE, parallel = parallel, 
                     pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(!missing(DoF) && missing(range) && missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>500){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          nllh_estfjw(par=param, z = z, lDoF = log(DoF), sites = sites, hit = hit, 
                      acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                      split = FALSE, parallel = parallel, 
                      pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(missing(DoF) && !missing(range) && missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          nllh_estrjw(par=param, z = z, range = range, sites = sites, hit = hit, 
                      acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                      split = FALSE, parallel = parallel, 
                      pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(missing(DoF) && missing(range) && !missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(exp(param[2])>500){return(Inf)}
          nllh_estsjw(par=param, z = z, smooth = smooth, sites = sites, hit = hit, 
                      acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                      split = FALSE, parallel = parallel, 
                      pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(missing(DoF) && missing(range) && missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(exp(param[2])>500){return(Inf)}
          if(round(inv.logit(param[3]),6)==0 || round(inv.logit(param[3]),6)==1){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_esta0jw(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                             acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_esta1jw(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                             acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_esta2jw(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                             acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_esta01jw(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                              acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_esta02jw(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                              acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_esta12jw(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                              acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_esta012jw(par=param, z = z, alpha = alpha, sites = sites, hit = hit, 
                           acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                           split = FALSE, parallel = parallel, 
                           pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
        }
      }else if(!missing(DoF) && !missing(range) && missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(round(inv.logit(param[1]),6)==0 || round(inv.logit(param[1]),6)==1){return(Inf)}
          nllh_estfrjw(par=param, z = z, lDoF = log(DoF), range = range, sites = sites, hit = hit, 
                       acov1=acov1, acov2=acov2, jw = jw, cmat = cmat,
                       split = FALSE, parallel = parallel, 
                       pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(!missing(DoF) && missing(range) && !missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>500){return(Inf)}
          nllh_estfsjw(par=param, z = z, lDoF = log(DoF), smooth = smooth, sites = sites, hit = hit, 
                       acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                       split = FALSE, parallel = parallel, 
                       pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(!missing(DoF) && missing(range) && missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estfa0jw(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estfa1jw(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estfa2jw(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estfa01jw(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estfa02jw(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estfa12jw(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estfa012jw(par=param, z = z, lDoF = log(DoF),alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
        }
      }else if(missing(DoF) && !missing(range) && !missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          nllh_estrsjw(par=param, z = z, range = range, smooth = smooth, sites = sites, hit = hit, 
                       acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                       split = FALSE, parallel = parallel, 
                       pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(missing(DoF) && !missing(range) && missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(round(inv.logit(param[2]),6)==0 || round(inv.logit(param[2]),6)==1){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estra0jw(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estra1jw(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estra2jw(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estra01jw(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estra02jw(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estra12jw(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estra012jw(par=param, z = z, range = range, alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
          
        }
      }else if(missing(DoF) && missing(range) && !missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param[1])>50){return(Inf)}
          if(exp(param[2])>500){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estsa0jw(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estsa1jw(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estsa2jw(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                              acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                              split = FALSE, parallel = parallel, 
                              pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estsa01jw(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estsa02jw(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estsa12jw(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estsa012jw(par=param, z = z, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                            acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                            split = FALSE, parallel = parallel, 
                            pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
          
        }
      }else if(!missing(DoF) && !missing(range) && !missing(smooth) && missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          nllh_estfrsjw(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, sites = sites, hit = hit, 
                        acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                        split = FALSE, parallel = parallel, 
                        pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
        }
      }else if(!missing(DoF) && !missing(range) && missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(round(inv.logit(param),6)==0 || round(inv.logit(param),6)==1){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estfra0jw(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estfra1jw(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estfra2jw(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estfra01jw(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                                acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estfra02jw(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                                acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estfra12jw(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                                acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estfra012jw(par=param, z = z, lDoF = log(DoF), range = range, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
        }
      }else if(!missing(DoF) && missing(range) && !missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param)>500){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estfsa0jw(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estfsa1jw(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estfsa2jw(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estfsa01jw(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estfsa02jw(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                                acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estfsa12jw(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                                acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estfsa012jw(par=param, z = z, lDoF = log(DoF), smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
        }
      }else if(missing(DoF) && !missing(range) && !missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
          if(exp(param)>50){return(Inf)}
          
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estrsa0jw(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estrsa1jw(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estrsa2jw(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                               acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                               split = FALSE, parallel = parallel, 
                               pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estrsa01jw(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                                acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estrsa02jw(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                                acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estrsa12jw(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                                acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          }else{ # all alpha0, alpha1, alpha2 are fixed
            nllh_estrsa012jw(par=param, z = z, range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                             acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                             split = FALSE, parallel = parallel, 
                             pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
          }
          
        }
      }else if(!missing(DoF) && !missing(range) && !missing(smooth) && !missing(alpha)){
        nllh_tmp <- function(param, z=z, hit=hit){
    
          if(any(is.na(alpha))){
            ind <- which(!is.na(alpha))
            alpha <- na.omit(alpha)
            if(length(ind)==1){
              if(ind==1){
                nllh_estfrsa0jw(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                                acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==2){
                nllh_estfrsa1jw(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                                acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind==3){
                nllh_estfrsa2jw(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                                acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                split = FALSE, parallel = parallel, 
                                pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              } 
            }else if(length(ind)==2){
              if(ind[1]==1 && ind[2]==2){
                nllh_estfrsa01jw(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                                 acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                 split = FALSE, parallel = parallel, 
                                 pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==1 && ind[2]==3){
                nllh_estfrsa02jw(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                                 acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                 split = FALSE, parallel = parallel, 
                                 pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }else if(ind[1]==2 && ind[2]==3){
                nllh_estfrsa12jw(par=param, z = z, lDoF = log(DoF), range = range, smooth = smooth, alpha = alpha, sites = sites, hit = hit, 
                                 acov1 = acov1, acov2 = acov2, jw = jw, cmat = cmat,
                                 split = FALSE, parallel = parallel, 
                                 pfun = mypmvsext, args1 = args1, args2 = args2, seed = seed)
              }
            }
          } # END if any(is.na(alpha))
          
        } # END nllh_tmp
      } 
      
    } # END else jw != Ns 
    
  } # END if model == "EST"
  
  xx <- optim(par0, nllh_tmp, z=z, hit=hit, method=method, control=control)
  
  if(model == "ET"){
    
    if(missing(DoF) && missing(smooth) && missing(range)){
      est <- c(exp(xx$par[1:2]), 2*inv.logit(xx$par[3]))
      names(est) <- c("DoF", "range", "smooth")
    }
    
    if(missing(DoF) && missing(smooth) && !missing(range)){
      est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]))
      names(est) <- c("DoF", "smooth")
    }
    
    if(missing(DoF) && !missing(smooth) && missing(range)){
      est <- exp(xx$par)
      names(est) <- c("DoF", "range")
    }
    
    if(!missing(DoF) && missing(smooth) && missing(range)){
      est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]))
      names(est) <- c("range", "smooth")
    }
    
    if(!missing(DoF) && !missing(smooth) && missing(range)){
      est <- exp(xx$par)
      names(est) <- "smooth"
    }
    
    if(!missing(DoF) && missing(smooth) && !missing(range)){
      est <- 2*inv.logit(xx$par)
      names(est) <- "range"
    }
    
    if(missing(DoF) && !missing(smooth) && !missing(range)){
      est <- exp(xx$par)
      names(est) <- "DoF"
    }
    
  } # END if model == "ET
  
  if(model == "EST"){
    
    if(missing(DoF) && missing(range) && missing(smooth) && missing(alpha) ){
      est <- c(exp(xx$par[1:2]), 2*inv.logit(xx$par[3]), xx$par[4:6])
      names(est) <- c("DoF", "range", "smooth", "alpha.0", "alpha.1", "alpha.2")
    }else if(!missing(DoF) && missing(range) && missing(smooth) && missing(alpha)){
      est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3:5])
      names(est) <- c("range", "smooth", "alpha.0", "alpha.1", "alpha.2")
    }else if(missing(DoF) && !missing(range) && missing(smooth) && missing(alpha)){
      est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3:5])
      names(est) <- c("DoF", "smooth", "alpha.0", "alpha.1", "alpha.2")
    }else if(missing(DoF) && missing(range) && !missing(smooth) && missing(alpha)){
      est <- c(exp(xx$par[1:2]), xx$par[3:5])
      names(est) <- c("DoF", "range", "alpha.0", "alpha.1", "alpha.2")
    }else if(missing(DoF) && missing(range) && missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(ind==1){
            est <- c(exp(xx$par[1:2]), 2*inv.logit(xx$par[3]), xx$par[4:5])
            names(est) <- c("DoF", "range", "smooth", "alpha1", "alpha2")
          }else if(ind==2){
            est <- c(exp(xx$par[1:2]), 2*inv.logit(xx$par[3]), xx$par[4:5])
            names(est) <- c("DoF", "range", "smooth", "alpha0", "alpha2")
          }else if(ind==3){
            est <- c(exp(xx$par[1:2]), 2*inv.logit(xx$par[3]), xx$par[4:5])
            names(est) <- c("DoF", "range", "smooth", "alpha1", "alpha1")
          } 
        }else if(length(ind)==2){
          if(ind[1]==1 && ind[2]==2){
            est <- c(exp(xx$par[1:2]), 2*inv.logit(xx$par[3]), xx$par[4])
            names(est) <- c("DoF", "range", "smooth", "alpha2")
          }else if(ind[1]==1 && ind[2]==3){
            est <- c(exp(xx$par[1:2]), 2*inv.logit(xx$par[3]), xx$par[4])
            names(est) <- c("DoF", "range", "smooth", "alpha1")
          }else if(ind[1]==2 && ind[2]==3){
            est <- c(exp(xx$par[1:2]), 2*inv.logit(xx$par[3]), xx$par[4])
            names(est) <- c("DoF", "range", "smooth", "alpha0")
          }
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        est <- c(exp(xx$par[1:2]), 2*inv.logit(xx$par[3]))
        names(est) <- c("DoF", "range", "smooth")
      }
    }else if(!missing(DoF) && !missing(range) && missing(smooth) && missing(alpha)){
      est <- c(2*inv.logit(xx$par[1]), xx$par[2:4])
      names(est) <- c("smooth", "alpha.0", "alpha.1", "alpha.2")
    }else if(!missing(DoF) && missing(range) && !missing(smooth) && missing(alpha)){
      est <- c(exp(xx$par[1]), xx$par[2:4])
      names(est) <- c("range", "alpha.0", "alpha.1", "alpha.2")
    }else if(!missing(DoF) && missing(range) && missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(ind==1){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3:4])
            names(est) <- c("range", "smooth", "alpha1", "alpha2")
          }else if(ind==2){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3:4])
            names(est) <- c("range", "smooth", "alpha0", "alpha2")
          }else if(ind==3){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3:4])
            names(est) <- c("range", "smooth", "alpha0", "alpha1")
          } 
        }else if(length(ind)==2){
          if(ind[1]==1 && ind[2]==2){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3])
            names(est) <- c("range", "smooth", "alpha2")
          }else if(ind[1]==1 && ind[2]==3){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3])
            names(est) <- c("range", "smooth", "alpha1")
          }else if(ind[1]==2 && ind[2]==3){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3])
            names(est) <- c("range", "smooth", "alpha0")
          }
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]))
        names(est) <- c("range", "smooth")
      }
    }else if(missing(DoF) && !missing(range) && !missing(smooth) && missing(alpha)){
      est <- c(exp(xx$par[1]), xx$par[2:4])
      names(est) <- c("DoF", "alpha.0", "alpha.1", "alpha.2")
    }else if(missing(DoF) && !missing(range) && missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(ind==1){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3:4])
            names(est) <- c("DoF", "smooth", "alpha1", "alpha2")
          }else if(ind==2){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3:4])
            names(est) <- c("DoF", "smooth", "alpha0", "alpha2")
          }else if(ind==3){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3:4])
            names(est) <- c("DoF", "smooth", "alpha0", "alpha1")
          } 
        }else if(length(ind)==2){
          if(ind[1]==1 && ind[2]==2){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3])
            names(est) <- c("DoF", "smooth", "alpha2")
          }else if(ind[1]==1 && ind[2]==3){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3])
            names(est) <- c("DoF", "smooth", "alpha1")
          }else if(ind[1]==2 && ind[2]==3){
            est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]), xx$par[3])
            names(est) <- c("DoF", "smooth", "alpha2")
          }
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        est <- c(exp(xx$par[1]), 2*inv.logit(xx$par[2]))
        names(est) <- c("DoF", "smooth")
      }
    }else if(missing(DoF) && missing(range) && !missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(ind==1){
            est <- c(exp(xx$par[1:2]), xx$par[3:4])
            names(est) <- c("DoF", "range", "alpha1", "alpha2")
          }else if(ind==2){
            est <- c(exp(xx$par[1:2]), xx$par[3:4])
            names(est) <- c("DoF", "range", "alpha0", "alpha2")
          }else if(ind==3){
            est <- c(exp(xx$par[1:2]), xx$par[3:4])
            names(est) <- c("DoF", "range", "alpha0", "alpha1")
          } 
        }else if(length(ind)==2){
          if(ind[1]==1 && ind[2]==2){
            est <- c(exp(xx$par[1:2]), xx$par[3])
            names(est) <- c("DoF", "range", "alpha2")
          }else if(ind[1]==1 && ind[2]==3){
            est <- c(exp(xx$par[1:2]), xx$par[3])
            names(est) <- c("DoF", "range", "alpha1")
          }else if(ind[1]==2 && ind[2]==3){
            est <- c(exp(xx$par[1:2]), xx$par[3])
            names(est) <- c("DoF", "range", "alpha0")
          }
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        est <- c(exp(xx$par[1:2]))
        names(est) <- c("DoF", "range")
      }
    }else if(!missing(DoF) && !missing(range) && !missing(smooth) && missing(alpha)){
      est <- xx$par[1:3]
      names(est) <- c("alpha.0", "alpha.1", "alpha.2")
    }else if(!missing(DoF) && !missing(range) && missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(ind==1){
            est <- c(2*inv.logit(xx$par[1]), xx$par[2:3])
            names(est) <- c("smooth", "alpha1", "alpha2")
          }else if(ind==2){
            est <- c(2*inv.logit(xx$par[1]), xx$par[2:3])
            names(est) <- c("smooth", "alpha0", "alpha2")
          }else if(ind==3){
            est <- c(2*inv.logit(xx$par[1]), xx$par[2:3])
            names(est) <- c("smooth", "alpha0", "alpha1")
          } 
        }else if(length(ind)==2){
          if(ind[1]==1 && ind[2]==2){
            est <- c(2*inv.logit(xx$par[1]), xx$par[2])
            names(est) <- c("smooth", "alpha2")
          }else if(ind[1]==1 && ind[2]==3){
            est <- c(2*inv.logit(xx$par[1]), xx$par[2])
            names(est) <- c("smooth", "alpha1")
          }else if(ind[1]==2 && ind[2]==3){
            est <- c(2*inv.logit(xx$par[1]), xx$par[2])
            names(est) <- c("smooth", "alpha0")
          }
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        est <- c(2*inv.logit(xx$par[1]) )
        names(est) <- c("smooth")
      }
    }else if(!missing(DoF) && missing(range) && !missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(ind==1){
            est <- c(exp(xx$par[1]),xx$par[2:3])
            names(est) <- c("range", "alpha1", "alpha2")
          }else if(ind==2){
            est <- c(exp(xx$par[1]),xx$par[2:3])
            names(est) <- c("range", "alpha0", "alpha2")
          }else if(ind==3){
            est <- c(exp(xx$par[1]),xx$par[2:3])
            names(est) <- c("range", "alpha0", "alpha1")
          } 
        }else if(length(ind)==2){
          if(ind[1]==1 && ind[2]==2){
            est <- c(exp(xx$par[1]),xx$par[2])
            names(est) <- c("range", "alpha2")
          }else if(ind[1]==1 && ind[2]==3){
            est <- c(exp(xx$par[1]),xx$par[2])
            names(est) <- c("range", "alpha1")
          }else if(ind[1]==2 && ind[2]==3){
            est <- c(exp(xx$par[1]),xx$par[2])
            names(est) <- c("range", "alpha0")
          }
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        est <- c(exp(xx$par[1]))
        names(est) <- c("range")
      }
    }else if(missing(DoF) && !missing(range) && !missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(ind==1){
            est <- c(exp(xx$par[1]),xx$par[2:3])
            names(est) <- c("DoF", "alpha1", "alpha2")
          }else if(ind==2){
            est <- c(exp(xx$par[1]),xx$par[2:3])
            names(est) <- c("DoF", "alpha0", "alpha2")
          }else if(ind==3){
            est <- c(exp(xx$par[1]),xx$par[2:3])
            names(est) <- c("DoF", "alpha0", "alpha1")
          } 
        }else if(length(ind)==2){
          if(ind[1]==1 && ind[2]==2){
            est <- c(exp(xx$par[1]),xx$par[2])
            names(est) <- c("DoF", "alpha2")
          }else if(ind[1]==1 && ind[2]==3){
            est <- c(exp(xx$par[1]),xx$par[2])
            names(est) <- c("DoF", "alpha1")
          }else if(ind[1]==2 && ind[2]==3){
            est <- c(exp(xx$par[1]),xx$par[2])
            names(est) <- c("DoF", "alpha0")
          }
        }
      }else{ # all alpha0, alpha1, alpha2 are fixed
        est <- c(exp(xx$par[1]))
        names(est) <- c("DoF")
      }
    }else if(!missing(DoF) && !missing(range) && !missing(smooth) && !missing(alpha)){
      if(any(is.na(alpha))){
        ind <- which(!is.na(alpha))
        if(length(ind)==1){
          if(ind==1){
            est <- xx$par
            names(est) <- c("alpha1", "alpha2")
          }else if(ind==2){
            est <- xx$par
            names(est) <- c("alpha0", "alpha2")
          }else if(ind==3){
            est <- xx$par
            names(est) <- c("alpha0", "alpha1")
          } 
        }else if(length(ind)==2){
          if(ind[1]==1 && ind[2]==2){
            est <- xx$par
            names(est) <- c("alpha2")
          }else if(ind[1]==1 && ind[2]==3){
            est <- xx$par
            names(est) <- c("alpha1")
          }else if(ind[1]==2 && ind[2]==3){
            est <- xx$par
            names(est) <- c("alpha0")
          }
        } # END if length(ind)==2
      } # END if any(is.na(alpha))
    }  
    
  }
  
  ## Compute standard errors using sandwich information matrix
 
  if(sandwich){
    
    epsVec <- rep(0.005, length(xx$par))
    score <- matrix(nrow=nrow(z), ncol=length(xx$par))
    
    fun <- function(par){ nllh_tmp(param=par, z=z, hit=hit)}
    bread <- fdhessian(fun=fun, par=xx$par, epsVec=epsVec)
    
    for(i in 1:nrow(z)){
        funi <- function(par){ nllh_tmp(param=par, z=matrix(z[i,], nrow=1), hit=matrix(hit[i,],nrow=1))}
        score[i,] <- fdjacobian(fun = funi, par = xx$par, split = FALSE, epsVec = epsVec)
    }
    
    meat <- var(score)
    Sand <- solve(bread) %*% meat %*% solve(bread) * nrow(z)
    
    sand <- sqrt(diag(Sand))
    TIC <- -2 * ( -xx$value - matrix.trace(meat %*% solve(bread) * nrow(z)) )
  }
  
  if(parallel){
    gc()
    closeAllConnections()
  }
 
  if(jw < Ns){
    
    if(sandwich){
      return(list(est=est, jw=jw, cmat=cmat, LL=-xx$value, stderr.sand=sand, TIC=TIC))
    }else{
      return(list(est=est, jw=jw, cmat=cmat, LL=-xx$value))
    }

  }else{
    
    if(sandwich){
      return(list(est=est, jw=jw, LL=-xx$value, stderr.sand=sand, TIC=TIC))
    }else{
      return(list(est=est, jw=jw, LL=-xx$value))
    }
    
  }
  
}


###############################################################################
###############################################################################
## Hidden functions
###############################################################################
###############################################################################

# Functions to do logit and inverse logit transformations
logit <- function(p){ return(log(p/(1-p)))}
# inv.logit <- function(x){ return(exp(x)/(1+exp(x)))}
inv.logit <- function(x){ return(1/(1+exp(-x)))}

# Make symmetric
makeSymmetric <- function(mat) 
{
  if(!isSymmetric(mat))
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  if(any(diag(mat) <= 0)) {
    print(diag(mat))
    stop("non-positive diagonal in makeSymmetric")
  }
  mat
}

# converts hitting scenario from matrix to lists of lists
# can specify a subset of sites for j-wise likelihood
convertHit <- function(mat, siteid = 1:ncol(mat), drop = FALSE) 
{
  if(!is.matrix(mat)) {
    stop("argument mat should be a matrix")
  }
  ns <- ncol(mat)
  if(!is.atomic(siteid) || length(siteid) <= 1)
    stop("siteid should be a vector of length two or more")
  if(!is.numeric(siteid) || !all(siteid %in% 1:ns) || any(duplicated(siteid)))
    stop("siteid has unknown or duplicated sites")
  
  convert1 <- function(vec) {
    out <- list()
    vec <- as.integer(factor(vec))
    for(i in 1:max(vec)) {
      out[[i]] <- which(vec == i)
    }
    out
  }
  res <- apply(mat[,siteid, drop = FALSE], 1, convert1)
  if(drop && length(res) == 1) res <- res[[1]]
  res
}

# central finite difference jacobian evaluation
fdjacobian <- function(fun, par, split, epsVec) 
{
  np <- length(par)
  if(np != length(epsVec)) stop("wrong length for epsVec")
  
  parmatp <- parmatm <- matrix(par, nrow=np, ncol=np)
  for(i in 1:np) {
    parmatp[i,i] <- parmatp[i,i] + epsVec[i]
    parmatm[i,i] <- parmatm[i,i] - epsVec[i]
  }
  
  #parp0 <- c(par[1] + eps1, par[2])
  #parm0 <- c(par[1] - eps1, par[2])
  #par0p <- c(par[1], par[2] + eps2)
  #par0m <- c(par[1], par[2] - eps2)
  #jmat <- matrix(NA, nrow = nrow(z), ncol = np)
  jmat <- matrix(NA,  ncol = np)
  for(i in 1:np)
    jmat[,i] <- (fun(parmatp[,i]) - 
                   fun(parmatm[,i]))/(2 * epsVec[i])
  if(!split) jmat <- colSums(jmat)
  jmat
}

# central finite difference hessian evaluation
fdhessian <- function(fun, par, epsVec) 
{
  np <- length(par)
  if(np != length(epsVec)) stop("wrong length for epsVec")
  
  
  fpar <- fun(par)
  parmatp <- parmatm <- matrix(par, nrow=np, ncol=np)
  for(i in 1:np) {
    parmatp[i,i] <- parmatp[i,i] + epsVec[i]
    parmatm[i,i] <- parmatm[i,i] - epsVec[i]
  }
  cp <- combn(np,2)
  parmatpp <- parmatmm <- parmatmp <- parmatpm <- matrix(par, nrow=np, ncol=ncol(cp))
  for(i in 1:ncol(cp)) {
    cp1 <- cp[1,i]; cp2 <- cp[2,i]
    parmatpp[cp1,i] <- parmatpp[cp1,i] + epsVec[cp1]
    parmatpp[cp2,i] <- parmatpp[cp2,i] + epsVec[cp2]
    parmatmm[cp1,i] <- parmatmm[cp1,i] - epsVec[cp1]
    parmatmm[cp2,i] <- parmatmm[cp2,i] - epsVec[cp2]
    parmatpm[cp1,i] <- parmatpm[cp1,i] + epsVec[cp1]
    parmatpm[cp2,i] <- parmatpm[cp2,i] - epsVec[cp2]
    parmatmp[cp1,i] <- parmatmp[cp1,i] - epsVec[cp1]
    parmatmp[cp2,i] <- parmatmp[cp2,i] + epsVec[cp2]
  }
  #parp0 <- c(par[1] + eps1, par[2])
  #parm0 <- c(par[1] - eps1, par[2])
  #par0p <- c(par[1], par[2] + eps2)
  #par0m <- c(par[1], par[2] - eps2)
  
  #parpp <- c(par[1] + eps1, par[2] + eps2)
  #parpm <- c(par[1] + eps1, par[2] - eps2)
  #parmp <- c(par[1] - eps1, par[2] + eps2)
  #parmm <- c(par[1] - eps1, par[2] - eps2)
  hmat <- matrix(NA, nrow = np, ncol = np)
  for(i in 1:np) 
    hmat[i,i] <- (fun(parmatp[,i]) + fun(parmatm[,i]) - 2*fpar)/(epsVec[i] * epsVec[i])
  for(i in 1:ncol(cp)) {
    cp1 <- cp[1,i]; cp2 <- cp[2,i]
    hmat[cp1,cp2] <- hmat[cp2,cp1] <- (fun(parmatpp[,i]) + fun(parmatmm[,i]) - 
                                         fun(parmatmp[,i]) - fun(parmatpm[,i])) / (4 * epsVec[cp1] * epsVec[cp2])
  }
  
  #h11 <- (fun(parp0, ...) + fun(parm0, ...) - 2*fpar)/(eps1 * eps1)
  #h22 <- (fun(par0p, ...) + fun(par0m, ...) - 2*fpar)/(eps2 * eps2)
  #h12 <- (fun(parpp, ...) + fun(parmm, ...) - fun(parmp, ...) - fun(parpm, ...)) / (4 * eps1 * eps2)
  #matrix(c(h11,h12,h12,h22), ncol = 2, nrow = 2)
  hmat
}

# fdsandwich <- function(fun, par, split, epsVec) 
# {
#   bread <- solve(fdhessian(fun = fun, par = par, epsVec = epsVec))
#   meat <- fdjacobian(fun = fun, par = par, split = split, epsVec = epsVec)
#   meat <- lapply(as.data.frame(t(meat)), function(x) x %o% x)
#   meat <- Reduce("+", meat)/length(meat)
#   bread %*% meat %*% bread
# }


getjwise <- function(sites, jw, thresh) 
{
  if(ncol(sites) != 2) 
    stop("sites must have two columns")
  ns <- nrow(sites)
  if(length(jw) != 1 || !is.numeric(jw))
    stop("jw must be a single value")
  jw <- as.integer(jw)
  if(!(jw %in% 2:ns))
    stop("jw must be between 2 and the number of sites")
  
  hmat <- as.matrix(dist(sites))
  out <- matrix(NA, nrow = jw, ncol = 0)
  
  for(i in 1:(ns-jw+1)) {
    within <- which(hmat[i,(i+1):ns] < thresh) + i
    if(length(within) == 1 && jw == 2) {
      out <- cbind(out, matrix(c(i, within), ncol = 1))
    } else if(length(within) >= (jw-1)) {
      out <- cbind(out, rbind(i, combn(within, jw-1), deparse.level = 0))
    }
  }
  
  out
}

###############################################################################
## NORMAL CDF
###############################################################################

mypmvnorm <- function(upper, sigma=diag(d), 
                      Nmax = 5000L, Nmin = 50L, eps = 0.0001, logeps = FALSE)
{
  d <- length(upper)
  if(nrow(sigma) != d || ncol(sigma) != d)
    stop("sigma is wrong dimension")
  cmat <- chol(sigma)
  
  .C("mypmvnorm", as.double(upper), as.integer(d), as.double(cmat),
     as.integer(Nmax), as.integer(Nmin), as.double(eps), as.integer(logeps),
     out = double(1))$out
}

###############################################################################
## STUDENT-T CDF
###############################################################################

mypmvt <- function(upper, sigma=diag(d), df = 1,
                   Nmax = 5000L, Nmin = 50L, eps = 0.0001, logeps = FALSE)
{
  d <- length(upper)
  if(nrow(sigma) != d || ncol(sigma) != d)
    stop("sigma is wrong dimension")
  cmat <- chol(sigma)
  
  .C("mypmvt", as.double(upper), as.integer(d), as.double(cmat), as.double(df),
     as.integer(Nmax), as.integer(Nmin), as.double(eps), as.integer(logeps),
     out = double(1))$out
}
###############################################################################
## skew extremal t distribution function
## does not have mean or non-centrality parameters
###############################################################################

mypmvsext <- function(upper, sigma=diag(d), df = 1, alpha = rep(0,d), ext = 0,
                      Nmax = 5000L, Nmin = 50L, eps = 0.0001, logeps = FALSE)
{
  d <- length(upper)
  qf <- function(x, A) colSums(x * (A %*% x))
  if(nrow(sigma) != d || ncol(sigma) != d)
    stop("sigma is wrong dimension")
  if(length(alpha) != d)
    stop("alpha is wrong length")
  if(length(ext) != 1)
    stop("ext must be a single value")
  
  upper <- upper / sqrt(diag(sigma)) 
  dimat <- diag(1/sqrt(diag(sigma)), nrow = nrow(sigma))
  rmat <- dimat %*% sigma %*% dimat
  den <- sqrt(1 + qf(alpha, rmat))
  extb <- ext / den; upper <- c(upper, extb)
  deltavec <- (rmat %*% alpha) / den
  rmatstar <- rbind(cbind(rmat, -deltavec), c(-deltavec, 1))
  cmat <- chol(rmatstar)
  
  out <- .C("mypmvt", as.double(upper), as.integer(d+1), as.double(cmat), as.double(df),
            as.integer(Nmax), as.integer(Nmin), as.double(eps), as.integer(logeps),
            out = double(1))$out
  out / pt(extb, df = df)
}


###############################################################################
## Routines for Extremal-t model 
###############################################################################

# correlation model is rho(h) = exp(-(h/lambda)^smooth)
# take par as log(DoF) and log(lambda)
# smooth is fixed
# DoF = 1 for Schalther (note degrees of freedon is actually DoF+1)
nllh_ets <- function(par, z, sites, hit, smooth = 1, split = FALSE, parallel = FALSE, 
                     pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{
  set.seed(seed)
  if(!is.numeric(z) || !is.matrix(z)) 
    stop("z must be a numeric matrix")
  nn <- nrow(z)
  if(any(z < 0)) stop("margins must be standard Frechet")
  
  hit <- convertHit(hit) 
  if(!is.list(hit) || length(hit) != nn)
    stop("incorrect format for hitting scenarios")
  
  par[1] <- exp(par[1]) # DoF
  par[2] <- exp(par[2]) # lambda
  
  if(par[1] < 0.1 || par[1] > 100) {
    return(Inf)
  }
  if(par[2] < 0.1 || par[2] > 100) {
    return(Inf)
  }
  
  if(parallel) {
    
    nllh <- foreach(i=1:nn, .combine = c, .export = c("vfun_et","vdfun_et","makeSymmetric")) %dopar% {
      v1 <- vfun_et(z = z[i,], sites = sites, DoF = par[1], lambda = par[2], smooth = smooth, pfun = pfun, args = args1, seed = seed)
      v2 <- vdfun_et(z = z[i,], sites = sites, DoF = par[1], lambda = par[2], smooth = smooth, slst = hit[[i]], pfun = pfun, args = args2, seed = seed)
      v1 - v2
    }
    
  } else {
    
    nllh <- numeric(nn)
    for(i in 1:nn) {
      v1 <- vfun_et(z = z[i,], sites = sites, DoF = par[1], lambda = par[2], smooth = smooth, pfun = pfun, args = args1, seed = seed)
      v2 <- vdfun_et(z = z[i,], sites = sites, DoF = par[1], lambda = par[2], smooth = smooth, slst = hit[[i]], pfun = pfun, args = args2, seed = seed)
      nllh[i] <- (v1 - v2)
    }
  }
  if(!split) nllh <- sum(nllh)
  
  nllh
}

# correlation model is rho(h) = exp(-(h/lambda)^smooth)
# take par as log(DoF), log(lambda), logit(smooth/2) where logit(x) = log(x / (1-x))
# DoF = 1 for Schlather (note degrees of freedom is actually DoF+1)
nllh_et <- function(par, z, sites, hit, split = FALSE, parallel = FALSE, 
                    pfun = mypmvt, args1 = list(), args2 = list(), seed = 123){
  
  smooth <- 2* inv.logit(par[3])
  par <- par[1:2]
  nllh_ets(par = par, z = z, sites = sites, hit = hit, smooth = smooth, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
  
} 

# Estimates smooth and range, fixes DoF
# par is log(range), logit(smooth/2) where logit(x) = log(x / (1-x))
# lDoF is log(DoF)
nllh_etf <- function(par, lDoF, z, sites, hit, split = FALSE, parallel = FALSE, 
                     pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{
  smooth <- 2* inv.logit(par[2])
  par <- c(lDoF, par[1])
  nllh_ets(par = par, z = z, sites = sites, hit = hit, smooth = smooth, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# Estimates smooth and DoF, fixes range
# par is log(DoF), logit(smooth/2) where logit(x) = log(x / (1-x))
# lDoF is log(DoF)
nllh_etr <- function(par, range, z, sites, hit, split = FALSE, parallel = FALSE, 
                     pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{
  smooth <- 2* inv.logit(par[2])
  par <- c(par[1], log(range))
  nllh_ets(par = par, z = z, sites = sites, hit = hit, smooth = smooth, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}


# Fixed smooth and DoF
# par is log(range)
# lDoF is log(DoF)
nllh_etsf <- function(par, lDoF, z, sites, hit, smooth, split = FALSE, parallel = FALSE, 
                      pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{
  par <- c(lDoF, par)
  nllh_ets(par = par, z = z, sites = sites, hit = hit, smooth = smooth, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# Fixed smooth and range
# par is log(DoF)
# lDoF is log(DoF)
nllh_etsr <- function(par, range, z, sites, hit, smooth, split = FALSE, parallel = FALSE, 
                      pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{
  par <- c(par, log(range))
  nllh_ets(par = par, z = z, sites = sites, hit = hit, smooth = smooth, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# Fixed range and DoF
# par is logit(smooth/2) where logit(x) = log(x / (1-x))
# lDoF is log(DoF)
nllh_etfr <- function(par, range, lDoF, z, sites, hit, split = FALSE, parallel = FALSE, 
                      pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{
  smooth <- 2* inv.logit(par)
  nllh_ets(par = c(lDoF, log(range)), z = z, sites = sites, hit = hit, smooth = smooth, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}


#j-wise likelihood
nllh_etsjw <- function(par, z, sites, hit, smooth = 1, jw = 2, cmat,
                      split = FALSE, parallel = FALSE, pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{
  njw <- ncol(cmat)
  ny <- nrow(z)
  if(nrow(cmat) != jw)
    stop("cmat has incorrect number of rows")
  
  if(!split) {
    
    if(parallel) {
      
      nllhjw <- foreach(i=1:njw, .combine = c, .packages='foreach', .export = c("nllh_et","convertHit","vfun_et","vdfun_et","makeSymmetric")) %dopar% {
        snum <- cmat[,i]
        nllh_ets(par, matrix(z[,snum], ncol=jw), sites[snum,], matrix(hit[,snum], ncol=jw), smooth = smooth, split = split, 
                 parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
        
      }
      
    } else {
      
      nllhjw <- rep(NA, njw)
      for(i in 1:njw) {
        snum <- cmat[,i]
        nllhjw[i] <- nllh_ets(par, matrix(z[,snum], ncol=jw), sites[snum,], matrix(hit[,snum], ncol=jw), smooth = smooth, split = split, 
                              parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
      }
    }
    nllhjw <- sum(nllhjw)
    
  } else {
    
    if(parallel) {
      
      nllhjw <- foreach(i=1:njw, .combine = cbind, .packages='foreach', .export = c("nllh_et","convertHit","vfun_et","vdfun_et","makeSymmetric")) %dopar% {
        snum <- cmat[,i]
        nllh_ets(par, z[,snum], sites[snum,], hit[,snum], smooth = smooth, split = split, 
                 parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
        
      }
      
    } else {
      
      nllhjw <- matrix(rep(NA, njw*ny), nrow = ny, ncol = njw) 
      for(i in 1:njw) {
        snum <- cmat[,i]
        nllhjw[,i] <- nllh_ets(par, z[,snum], sites[snum,], hit[,snum], smooth = smooth, split = split, 
                               parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
      }
    }
    nllhjw <- rowSums(nllhjw)
    
  }
  
  nllhjw
}

# correlation model is rho(h) = exp(-(h/lambda)^smooth)
# take par as log(DoF), log(lambda), logit(smooth/2) where logit(x) = log(x / (1-x))
# DoF = 1 for Schalther (note degrees of freedon is actually DoF+1)
nllh_etjw <- function(par, z, sites, hit, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                      pfun = mypmvt, args1 = list(), args2 = list(), seed = 123){
  
  smooth <- 2* inv.logit(par[3])
  par <- par[1:2]
  nllh_etsjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
  
} 

nllh_etfjw <- function(par, lDoF, z, sites, hit, jw = 2, cmat,
                       split = FALSE, parallel = FALSE, pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{  
  smooth <- 2* inv.logit(par[2])
  par <- c(lDoF, par[1])
  nllh_etsjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}  

nllh_etrjw <- function(par, range, z, sites, hit, jw = 2, cmat,
                       split = FALSE, parallel = FALSE, pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{  
  smooth <- 2* inv.logit(par[2])
  par <- c(par[1], log(range))
  nllh_etsjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}  


nllh_etsfjw <- function(par, lDoF, z, sites, hit, smooth, jw = 2, cmat,
                        split = FALSE, parallel = FALSE, pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{  
  par <- c(lDoF, par)
  nllh_etsjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}  

nllh_etsrjw <- function(par, range, z, sites, hit, smooth, jw = 2, cmat,
                        split = FALSE, parallel = FALSE, pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{  
  par <- c(par, log(range))
  nllh_etsjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}  

nllh_etfrjw <- function(par, lDoF, range, z, sites, hit, jw = 2, cmat,
                        split = FALSE, parallel = FALSE, pfun = mypmvt, args1 = list(), args2 = list(), seed = 123) 
{  
  smooth <- 2* inv.logit(par)
  nllh_etsjw(par = c(lDoF, log(range)), z = z, sites = sites, hit = hit, smooth = smooth, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}  



# calculates V(z)
vfun_et <- function(z, sites, DoF, lambda, smooth,
                    split = FALSE, pfun = mypmvt, args = list(), seed) 
{
  set.seed(seed)
  d <- length(z)
  
  if(ncol(sites) != 2) stop("sites must have two columns")
  if(nrow(sites) != d) stop("incorrect number of sites")
  hmat <- as.matrix(dist(sites))
  
  # correlation model is rho(h) = exp(-(h/lambda)^smooth)
  rhomat <- exp(-(hmat/lambda)^smooth)
  rhomat2 <- sqrt(1 - rhomat*rhomat)
  
  vec <- numeric(d)
  for(i in 1:d) {
    rhoveci <- rhomat[i,-i]
    rhoveci2 <- rhomat2[i,-i]
    rhomati <- rhomat[-i,-i]
    
    eta <- sqrt(DoF+1) * ((z[i]/z[-i])^(-1/DoF) - rhoveci) / rhoveci2 
    rmat <- (rhomati - (rhoveci %o% rhoveci)) / (rhoveci2 %o% rhoveci2)
    Eigen <- tryCatch(eigen(rmat)$values, error=function(e) -1) 
    if(any(Eigen<=0)){return(-1e300)}
    
    vec[i] <- do.call(pfun, c(list(upper=eta, sigma=rmat, df = DoF+1), args))
  }
  
  if(split) return(vec/z)
  return(sum(vec/z))
}

# calculates sum of log(-V_s(z)) over s in hitting scenario
# slst is a list containing the hitting scenario
vdfun_et <- function(z, sites, DoF, lambda, smooth, slst = list(1:d), 
                     log.p = TRUE, split = FALSE, pfun = mypmvt, args = list(), seed) 
{
  qf <- function(x, A) colSums(x * (A %*% x))
  ldetinv <- function(mat) {
    cfac <- chol(mat)
    list(ldet = 2*sum(log(diag(cfac))), inv = chol2inv(cfac))
  }
  
  set.seed(seed)
  d <- length(z)
  if(ncol(sites) != 2) stop("sites must have two columns")
  if(nrow(sites) != d) stop("incorrect number of sites")
  tst <- sort(unlist(slst))
  if(length(tst) != d || !all(tst == 1:d)) 
    stop("invalid hitting scenario")
  
  hmat <- as.matrix(dist(sites))
  
  # correlation model is rho(h) = exp(-(h/lambda)^smooth)
  rhomat <- exp(-(hmat/lambda)^smooth)
  
  # log(cv)
  logcv <- log(pi)/2 - (DoF - 2)/2 * log(2) - lgamma((DoF + 1)/2)
  
  # variances are one
  sigma <- rhomat 
  isigma <- solve(sigma)
  
  # calculates log(-V_{sv}) for log.p TRUE
  vsfun <- function(z, sites, sv = 1:d, log.p = TRUE) 
  {
    s <- length(sv)
    zs <- z[sv]; zo <- z[-sv]
    
    sigmas <- sigma[sv, sv, drop = FALSE]
    ldi <- ldetinv(sigmas)
    ldsigmas <- ldi$ldet; isigmas <- ldi$inv
    afun <- qf(zs^(1/DoF), isigmas)
    
    sigmaso <- sigma[-sv, -sv, drop = FALSE]
    sigmaxx <- sigma[-sv, sv, drop = FALSE]
    
    # if sv is 1:d then mypmvt is taken as unity
    if(d != s) {
      tmp <- sigmaxx %*% isigmas 
      
      Sig <- sigmaso - tmp %*% t(sigmaxx)
      if(any(diag(Sig)<=0)){return(-1e300)}
      Gam <- afun / (s + DoF) * makeSymmetric(Sig)
      # Gam <- afun / (s + DoF) * makeSymmetric(sigmaso - tmp %*% t(sigmaxx))
      Eigen <- tryCatch(eigen(Gam)$values, error=function(e) -1) 
      if(any(Eigen<=0)){return(-1e300)}
      
      mu <- as.numeric(tmp %*% zs^(1/DoF))
      #p1 <- as.numeric(log(do.call(pfun, c(list(upper=zo^(1/DoF) - mu, sigma=Gam, df = DoF+s), args))))
      # to avoid bug in mvtnorm (reported so should be fixed soon)
      if(length(zo)==1){
        p1 <- as.numeric(log(do.call(pfun, c(list(upper=as.numeric((zo^(1/DoF) - mu)/sqrt(Gam)), df = DoF+s), args))))
      }else{
        p1 <- as.numeric(log(do.call(pfun, c(list(upper=zo^(1/DoF) - mu, sigma=Gam, df = DoF+s), args)))) 
      }
    } else {
      p1 <- 0
    }
    
    p2 <- logcv - (s-1) * log(DoF) + (DoF - 2)/2 * log(2) - s/2 * log(pi) - ldsigmas/2
    p3 <- lgamma((s + DoF)/2) + (1 - DoF)/DoF * sum(log(zs)) - (s + DoF)/2 * log(afun)   
    
    out <- p1 + p2 + p3
    if(!log.p) out <- exp(out)
    out
  }  
  
  val <- numeric(length(slst))
  for (k in 1:length(slst)) {
    val[k] <- vsfun(z = z, sites = sites, sv = slst[[k]], log.p = TRUE) 
  }
  if(!split) val <- sum(val)
  if(!log.p) val <- exp(val)
  val
}

###############################################################################
## Routines for Extremal Skew-t model 
###############################################################################

# correlation model is rho(h) = exp(-(h/lambda)^smooth)
# take par as log(DoF) and log(lambda)
# DoF = 1 for Schalther (note degrees of freedom is actually DoF+1)
# alpha is a skewness vector, default is zero.
nllh_EST <- function(par, z, sites, hit, smooth = 1, alpha = rep(0, nrow(sites)), split = FALSE, parallel = FALSE, 
                     pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  # set.seed(seed)
  if(!is.numeric(z) || !is.matrix(z)) 
    stop("z must be a numeric matrix")
  nn <- nrow(z)
  if(any(z < 0)) stop("margins must be standard Frechet")
  if(length(alpha) != nrow(sites)) stop("alpha is incorrect length")
  
  hit <- convertHit(hit) 
  if(!is.list(hit) || length(hit) != nn)
    stop("incorrect format for hitting scenarios")
  
  par[1] <- exp(par[1]) # DoF
  par[2] <- exp(par[2]) # lambda
  
  
  if(par[1] < 0.1 || par[1] > 100) {
    #browser()
    return(Inf)
  }
  if(par[2] < 0.0000001 || par[2] > 100) {
    #browser()
    return(Inf)
  }
  
  if(parallel) {
    
    nllh <- foreach(i=1:nn, .combine = c, .export = c("vfun_est","vdfun_est","makeSymmetric")) %dopar% {
      v1 <- vfun_est(z = z[i,], sites = sites, DoF = par[1], lambda = par[2], smooth = smooth, alpha = alpha, pfun = pfun, args = args1, seed = seed)
      v2 <- vdfun_est(z = z[i,], sites = sites, DoF = par[1], lambda = par[2], smooth = smooth, alpha = alpha, slst = hit[[i]], pfun = pfun, args = args2, seed = seed)
      v1 - v2
    }
    
  } else {
    
    nllh <- numeric(nn)
    for(i in 1:nn) {
      v1 <- vfun_est(z = z[i,], sites = sites, DoF = par[1], lambda = par[2], smooth = smooth, alpha = alpha, pfun = pfun, args = args1, seed = seed)
      v2 <- vdfun_est(z = z[i,], sites = sites, DoF = par[1], lambda = par[2], smooth = smooth, alpha = alpha, slst = hit[[i]], pfun = pfun, args = args2, seed = seed)
      nllh[i] <- (v1 - v2)
    }
  }
  if(!split) nllh <- sum(nllh)
  
  nllh
}

# par is (log(DoF), log(range), logit(smooth/2), skintercept, sk1, sk2)
nllh_est <- function(par, z, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                     pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{

  smooth <- 2 * inv.logit(par[3])
  Alpha <- par[4] + par[5]*acov1 + par[6]*acov2
  par <- par[1:2]
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# par is (log(range), logit(smooth/2), skintercept, sk1, sk2)
nllh_estf <- function(par, z, lDoF, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                     pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + par[4]*acov1 + par[5]*acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# par is (log(DoF), logit(smooth/2), skintercept, sk1, sk2)
nllh_estr <- function(par, z, range, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                      pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + par[4]*acov1 + par[5]*acov2
  par <- c(par[1], range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# par is (log(DoF), log(range), skintercept, sk1, sk2)
nllh_ests <- function(par, z, smooth, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                      pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[3] + par[4]*acov1 + par[5]*acov2
  par <- par[1:2]
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha0 is fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha1, alpha2)
nllh_esta0 <- function(par, z, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                       pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- alpha + par[4] * acov1 + par[5] * acov2
  par <- par[1:2]
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha1 is fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha0, alpha2)
nllh_esta1 <- function(par, z, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                       pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- par[4] + alpha * acov1 + par[5] * acov2
  par <- par[1:2]
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha2 is fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha0, alpha1)
nllh_esta2 <- function(par, z, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                       pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- par[4] + par[5] * acov1 + alpha * acov2
  par <- par[1:2]
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha0 and alpha1 are fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha2)
nllh_esta01 <- function(par, z, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                       pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- alpha[1] + alpha[2] * acov2 + par[4] * acov2
  par <- par[1:2]
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha0 and alpha2 are fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha1)
nllh_esta02 <- function(par, z, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- alpha[1] + par[4] * acov1 + alpha[2] * acov2
  par <- par[1:2]
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha1 and alpha2 are fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha0)
nllh_esta12 <- function(par, z, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- par[4] + alpha[1] * acov1 + alpha[2] * acov2
  par <- par[1:2]
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha0, alpha1 and alpha2 are fixed
# par is (log(DoF), log(range), logit(smooth/2))
nllh_esta012 <- function(par, z, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  par <- par[1:2]
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}


# lDoF = log(DoF) is fixed
# range is fixed
# par is (logit(smooth/2), skintercept, sk1, sk2)
nllh_estfr <- function(par, z, lDoF, range, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                      pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  alpha <- par[2] + par[3]*acov1 + par[4]*acov2
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# par is (log(range), skintercept, sk1, sk2)
nllh_estfs <- function(par, z, lDoF, smooth, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                       pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  alpha <- par[2] + par[3]*acov1 + par[4]*acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# alpha0 is fixed
# par is (log(range), logit(smooth/2), alpha1, alpha2 )
nllh_estfa0 <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                       pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha + par[3] * acov1 + par[4] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# alpha1 is fixed
# par is (log(range), logit(smooth/2), alpha0, alpha2 )
nllh_estfa1 <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + alpha * acov1 + par[4] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# alpha2 is fixed
# par is (log(range), logit(smooth/2), alpha0, alpha1 )
nllh_estfa2 <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + par[4] * acov1 + alpha * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# alpha0 and alpha1 are fixed
# par is (log(range), logit(smooth/2), alpha2 )
nllh_estfa01 <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + alpha[2] * acov1 + par[3] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# alpha0 and alpha2 are fixed
# par is (log(range), logit(smooth/2), alpha1 )
nllh_estfa02 <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + par[3] * acov1 + alpha[2] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# alpha1 and alpha2 are fixed
# par is (log(range), logit(smooth/2), alpha0 )
nllh_estfa12 <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + alpha[1] * acov1 + alpha[2] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (log(range), logit(smooth/2) )
nllh_estfa012 <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# par is (log(DoF), skintercept, sk1, sk2)
nllh_estrs <- function(par, z, range, smooth, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                       pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  alpha <- par[2] + par[3]*acov1 + par[4]*acov2
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha0 is fixed
# par is (log(DoF), logit(smooth/2), alpha1, alpha2 )
nllh_estra0 <- function(par, z, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha + par[3] * acov1 + par[4] * acov2
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha1 is fixed
# par is (log(DoF), logit(smooth/2), alpha0, alpha2 )
nllh_estra1 <- function(par, z, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + alpha * acov1 + par[4] * acov2
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha2 is fixed
# par is (log(DoF), logit(smooth/2), alpha0, alpha1 )
nllh_estra2 <- function(par, z, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + par[4] * acov1 + alpha * acov2
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha0 and alpha1 are fixed
# par is (log(DoF), logit(smooth/2), alpha2 )
nllh_estra01 <- function(par, z, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + alpha[2] * acov1 + par[3] * acov2
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha0 and alpha2 are fixed
# par is (log(DoF), logit(smooth/2), alpha1 )
nllh_estra02 <- function(par, z, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + par[3] * acov1 + alpha[2] * acov2
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha1 and alpha2 are fixed
# par is (log(DoF), logit(smooth/2), alpha0 )
nllh_estra12 <- function(par, z, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + alpha[1] * acov1 + alpha[2] * acov2
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (log(DoF), logit(smooth/2) )
nllh_estra012 <- function(par, z, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha0 is fixed
# par is (log(DoF), log(range), alpha1, alpha2)
nllh_estsa0 <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- alpha + par[3] * acov1 + par[4] * acov2
  nllh_EST(par = par[1:2], z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha1 is fixed
# par is (log(DoF), log(range), alpha0, alpha2)
nllh_estsa1 <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- par[3] + alpha * acov1 + par[4] * acov2
  nllh_EST(par = par[1:2], z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha2 is fixed
# par is (log(DoF), log(range), alpha0, alpha1)
nllh_estsa2 <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- par[3] + par[4] * acov1 + alpha * acov2
  nllh_EST(par = par[1:2], z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha0 and alpha1 are fixed
# par is (log(DoF), log(range), alpha2)
nllh_estsa01 <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- alpha[1] + alpha[2] * acov1 + par[3] * acov2
  nllh_EST(par = par[1:2], z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha0 and alpha2 are fixed
# par is (log(DoF), log(range), alpha1)
nllh_estsa02 <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- alpha[1] + par[3] * acov1 + alpha[2] * acov2
  nllh_EST(par = par[1:2], z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha1 and alpha2 are fixed
# par is (log(DoF), log(range), alpha0)
nllh_estsa12 <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- par[3] + alpha[1] * acov1 + alpha[2] * acov2
  nllh_EST(par = par[1:2], z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (log(DoF), log(range))
nllh_estsa012 <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  nllh_EST(par = par[1:2], z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# par is (skintercept, sk1, sk2)
nllh_estfrs <- function(par, z, lDoF, range, smooth, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                       pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  alpha <- par[1] + par[2]*acov1 + par[3]*acov2
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha0 is fixed
# par is (logit(smooth/2), alpha1, alpha2)
nllh_estfra0 <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- alpha + par[2] * acov1 + par[3]*acov2
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha1 is fixed
# par is (logit(smooth/2), alpha0, alpha2)
nllh_estfra1 <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- par[2] + alpha * acov1 + par[3]*acov2
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha2 is fixed
# par is (logit(smooth/2), alpha0, alpha1)
nllh_estfra2 <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- par[2] + par[3] * acov1 + alpha *acov2
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha0 and alpha1 are fixed
# par is (logit(smooth/2), alpha2)
nllh_estfra01 <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- alpha[1] + alpha[2] * acov1 + par[2] *acov2
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha0 and alpha2 are fixed
# par is (logit(smooth/2), alpha1)
nllh_estfra02 <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- alpha[1] + par[2] * acov1 + alpha[2] * acov2
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha1 and alpha2 are fixed
# par is (logit(smooth/2), alpha0)
nllh_estfra12 <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- par[2] + alpha[1] * acov1 + alpha[2] * acov2
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (logit(smooth/2))
nllh_estfra012 <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha0 is fixed
# par is (log(range), alpha1, alpha2)
nllh_estfsa0 <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha + par[2] * acov1 + par[3] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha1 is fixed
# par is (log(range), alpha0, alpha2)
nllh_estfsa1 <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + alpha * acov1 + par[3] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha2 is fixed
# par is (log(range), alpha0, alpha1)
nllh_estfsa2 <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + par[3] * acov1 + alpha * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha0 and alpha1 are fixed
# par is (log(range), alpha2)
nllh_estfsa01 <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + alpha[2] * acov1 + par[2] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha0 and alpha2 are fixed
# par is (log(range), alpha1)
nllh_estfsa02 <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + par[2] * acov1 + alpha[2] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha1 and alpha2 are fixed
# par is (log(range), alpha0)
nllh_estfsa12 <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + alpha[1] * acov1 + alpha[2] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (log(range))
nllh_estfsa012 <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  par <- c(lDoF,par[1])
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# alpha0 is fixed
# par is (log(DoF), alpha1, alpha2)
nllh_estrsa0 <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha + par[2] * acov1 + par[3] * acov2 
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# alpha1 is fixed
# par is (log(DoF), alpha0, alpha2)
nllh_estrsa1 <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + alpha * acov1 + par[3] * acov2 
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# alpha2 is fixed
# par is (log(DoF), alpha0, alpha1)
nllh_estrsa2 <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + par[3] * acov1 + alpha * acov2 
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# alpha0 and alpha1 are fixed
# par is (log(DoF), alpha2)
nllh_estrsa01 <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + alpha[2] * acov1 + par[2] * acov2 
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# alpha0 and alpha2 are fixed
# par is (log(DoF), alpha1)
nllh_estrsa02 <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + par[2] * acov1 + alpha[2] * acov2 
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# alpha1 and alpha2 are fixed
# par is (log(DoF), alpha0)
nllh_estrsa12 <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + alpha[1] * acov1 + alpha[2] * acov2 
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (log(DoF), alpha0)
nllh_estrsa012 <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2 
  par <- c(par[1],range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha0 is fixed
# par is (alpha1, alpha2)
nllh_estfrsa0 <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha + par[1] * acov1 + par[2] * acov2 
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha1 is fixed
# par is (alpha0, alpha2)
nllh_estfrsa1 <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[1] + alpha * acov1 + par[2] * acov2 
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha2 is fixed
# par is (alpha0, alpha1)
nllh_estfrsa2 <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[1] + par[2] * acov1 + alpha * acov2 
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha0 and alpha1 are fixed
# par is (alpha2)
nllh_estfrsa01 <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + alpha[2] * acov1 + par * acov2 
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha0 and alpha2 are fixed
# par is (alpha1)
nllh_estfrsa02 <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + par * acov1 + alpha[2] * acov2 
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha1 and alpha2 are fixed
# par is (alpha0)
nllh_estfrsa12 <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par + alpha[1] * acov1 + alpha[2] * acov2 
  par <- c(lDoF,range)
  nllh_EST(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = alpha, split = split, parallel = parallel, 
           pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

#j-wise likelihood
nllh_ESTjw <- function(par, z, sites, hit, smooth = 1, alpha = rep(0, nrow(sites)), jw = 2, cmat,
                       split = FALSE, parallel = FALSE, pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  njw <- ncol(cmat)
  ny <- nrow(z)
  if(nrow(cmat) != jw)
    stop("cmat has incorrect number of rows")
  
  if(!split) {
    
    if(parallel) {
      
      nllhjw <- foreach(i=1:njw, .combine = c, .packages='foreach', .export = c("nllh_est","convertHit","vfun_est","vdfun_est","makeSymmetric")) %dopar% {
        snum <- cmat[,i]
        nllh_EST(par, matrix(z[,snum], ncol=jw), sites[snum,], matrix(hit[,snum], ncol=jw), smooth = smooth, alpha = alpha[snum], split = FALSE, 
                 parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
        
      }
      
    } else {
      
      nllhjw <- rep(NA, njw)
      for(i in 1:njw) {
        snum <- cmat[,i]
        nllhjw[i] <- nllh_EST(par, matrix(z[,snum],ncol=jw), sites[snum,], matrix(hit[,snum],ncol=jw), smooth = smooth, alpha = alpha[snum], split = FALSE, 
                              parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
      }
    }
    nllhjw <- sum(nllhjw)
    
  } else {
    
    if(parallel) {
      
      nllhjw <- foreach(i=1:njw, .combine = cbind, .packages='foreach', .export = c("nllh_est","convertHit","vfun_est","vdfun_est", "makeSymmetric")) %dopar% {
        snum <- cmat[,i]
        nllh_EST(par, z[,snum], sites[snum,], hit[,snum], smooth = smooth, alpha = alpha[snum], split = TRUE, 
                 parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
        
      }
      
    } else {
      
      nllhjw <- matrix(rep(NA, njw*ny), nrow = ny, ncol = njw) 
      for(i in 1:njw) {
        snum <- cmat[,i]
        nllhjw[,i] <- nllh_EST(par, z[,snum], sites[snum,], hit[,snum], smooth = smooth, alpha = alpha[snum], split = TRUE, 
                               parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
      }
    }
    nllhjw <- rowSums(nllhjw)
    
  }
  
  nllhjw
}

# par is (log(DoF), log(range), logit(smooth/2), skintercept, sk1, sk2)
nllh_estjw <- function(par, z, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                       pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- par[4] + par[5]*acov1 + par[6]*acov2
  par <- par[1:2]
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# par is (log(range), logit(smooth/2), skintercept, sk1, sk2)
nllh_estfjw <- function(par, z, lDoF, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + par[4]*acov1 + par[5]*acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# par is (log(DoF), logit(smooth/2), skintercept, sk1, sk2)
nllh_estrjw <- function(par, z, range, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + par[4]*acov1 + par[5]*acov2
  par <- c(par[1], range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# par is (log(DoF), log(range), skintercept, sk1, sk2)
nllh_estsjw <- function(par, z, smooth, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[3] + par[4]*acov1 + par[5]*acov2
  par <- par[1:2]
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha0 is fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha1, alpha2)
nllh_esta0jw <- function(par, z, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- alpha + par[4] * acov1 + par[5] * acov2
  par <- par[1:2]
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)  
}

# alpha1 is fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha0, alpha2)
nllh_esta1jw <- function(par, z, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- par[4] + alpha * acov1 + par[5] * acov2
  par <- par[1:2]
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha2 is fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha0, alpha1)
nllh_esta2jw <- function(par, z, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- par[4] + par[5] * acov1 + alpha * acov2
  par <- par[1:2]
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha0 and alpha1 are fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha2)
nllh_esta01jw <- function(par, z, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- alpha[1] + alpha[2] * acov2 + par[4] * acov2
  par <- par[1:2]
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha0 and alpha2 are fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha1)
nllh_esta02jw <- function(par, z, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- alpha[1] + par[4] * acov1 + alpha[2] * acov2
  par <- par[1:2]
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha1 and alpha2 are fixed
# par is (log(DoF), log(range), logit(smooth/2), alpha0)
nllh_esta12jw <- function(par, z, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                        pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- par[4] + alpha[1] * acov1 + alpha[2] * acov2
  par <- par[1:2]
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# alpha0, alpha1 and alpha2 are fixed
# par is (log(DoF), log(range), logit(smooth/2))
nllh_esta012jw <- function(par, z, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[3])
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  par <- par[1:2]
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}


# lDoF = log(DoF) is fixed
# range is fixed
# par is (logit(smooth/2), skintercept, sk1, sk2)
nllh_estfrjw <- function(par, z, lDoF, range, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- par[2] + par[3]*acov1 + par[4]*acov2
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# par is (log(range), skintercept, sk1, sk2)
nllh_estfsjw <- function(par, z, lDoF, smooth, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + par[3]*acov1 + par[4]*acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# alpha0 is fixed
# par is (log(range), logit(smooth/2), alpha1, alpha2 )
nllh_estfa0jw <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha + par[3] * acov1 + par[4] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# alpha1 is fixed
# par is (log(range), logit(smooth/2), alpha0, alpha2 )
nllh_estfa1jw <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + alpha * acov1 + par[4] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# alpha2 is fixed
# par is (log(range), logit(smooth/2), alpha0, alpha1 )
nllh_estfa2jw <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + par[4] * acov1 + alpha * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# alpha0 and alpha1 are fixed
# par is (log(range), logit(smooth/2), alpha2 )
nllh_estfa01jw <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + alpha[2] * acov1 + par[3] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# alpha0 and alpha2 are fixed
# par is (log(range), logit(smooth/2), alpha1 )
nllh_estfa02jw <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + par[3] * acov1 + alpha[2] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# alpha1 and alpha2 are fixed
# par is (log(range), logit(smooth/2), alpha0 )
nllh_estfa12jw <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + alpha[1] * acov1 + alpha[2] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (log(range), logit(smooth/2) )
nllh_estfa012jw <- function(par, z, lDoF, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# par is (log(DoF), skintercept, sk1, sk2)
nllh_estrsjw <- function(par, z, range, smooth, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                         pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + par[3]*acov1 + par[4]*acov2
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha0 is fixed
# par is (log(DoF), logit(smooth/2), alpha1, alpha2 )
nllh_estra0jw <- function(par, z, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha + par[3] * acov1 + par[4] * acov2
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha1 is fixed
# par is (log(DoF), logit(smooth/2), alpha0, alpha2 )
nllh_estra1jw <- function(par, z, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + alpha * acov1 + par[4] * acov2
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# range is fixed
# alpha2 is fixed
# par is (log(DoF), logit(smooth/2), alpha0, alpha1 )
nllh_estra2jw <- function(par, z, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + par[4] * acov1 + alpha * acov2
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha0 and alpha1 are fixed
# par is (log(DoF), logit(smooth/2), alpha2 )
nllh_estra01jw <- function(par, z, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + alpha[2] * acov1 + par[3] * acov2
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha0 and alpha2 are fixed
# par is (log(DoF), logit(smooth/2), alpha1 )
nllh_estra02jw <- function(par, z, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + par[3] * acov1 + alpha[2] * acov2
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha1 and alpha2 are fixed
# par is (log(DoF), logit(smooth/2), alpha0 )
nllh_estra12jw <- function(par, z, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- par[3] + alpha[1] * acov1 + alpha[2] * acov2
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (log(DoF), logit(smooth/2) )
nllh_estra012jw <- function(par, z, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[2])
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha0 is fixed
# par is (log(DoF), log(range), alpha1, alpha2)
nllh_estsa0jw <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- alpha + par[3] * acov1 + par[4] * acov2
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha1 is fixed
# par is (log(DoF), log(range), alpha0, alpha2)
nllh_estsa1jw <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- par[3] + alpha * acov1 + par[4] * acov2
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha2 is fixed
# par is (log(DoF), log(range), alpha0, alpha1)
nllh_estsa2jw <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- par[3] + par[4] * acov1 + alpha * acov2
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# smooth is fixed
# alpha0 and alpha1 are fixed
# par is (log(DoF), log(range), alpha2)
nllh_estsa01jw <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- alpha[1] + alpha[2] * acov1 + par[3] * acov2
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha0 and alpha2 are fixed
# par is (log(DoF), log(range), alpha1)
nllh_estsa02jw <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- alpha[1] + par[3] * acov1 + alpha[2] * acov2
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha1 and alpha2 are fixed
# par is (log(DoF), log(range), alpha0)
nllh_estsa12jw <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- par[3] + alpha[1] * acov1 + alpha[2] * acov2
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# smooth is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (log(DoF), log(range))
nllh_estsa012jw <- function(par, z, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# par is (skintercept, sk1, sk2)
nllh_estfrsjw <- function(par, z, lDoF, range, smooth, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                          pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[1] + par[2]*acov1 + par[3]*acov2
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha0 is fixed
# par is (logit(smooth/2), alpha1, alpha2)
nllh_estfra0jw <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- alpha + par[2] * acov1 + par[3]*acov2
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha1 is fixed
# par is (logit(smooth/2), alpha0, alpha2)
nllh_estfra1jw <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- par[2] + alpha * acov1 + par[3]*acov2
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha2 is fixed
# par is (logit(smooth/2), alpha0, alpha1)
nllh_estfra2jw <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- par[2] + par[3] * acov1 + alpha *acov2
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha0 and alpha1 are fixed
# par is (logit(smooth/2), alpha2)
nllh_estfra01jw <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- alpha[1] + alpha[2] * acov1 + par[2] *acov2
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha0 and alpha2 are fixed
# par is (logit(smooth/2), alpha1)
nllh_estfra02jw <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- alpha[1] + par[2] * acov1 + alpha[2] * acov2
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha1 and alpha2 are fixed
# par is (logit(smooth/2), alpha0)
nllh_estfra12jw <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- par[2] + alpha[1] * acov1 + alpha[2] * acov2
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# range is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (logit(smooth/2))
nllh_estfra012jw <- function(par, z, lDoF, range, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                             pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  smooth <- 2 * inv.logit(par[1])
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha0 is fixed
# par is (log(range), alpha1, alpha2)
nllh_estfsa0jw <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha + par[2] * acov1 + par[3] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha1 is fixed
# par is (log(range), alpha0, alpha2)
nllh_estfsa1jw <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + alpha * acov1 + par[3] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha2 is fixed
# par is (log(range), alpha0, alpha1)
nllh_estfsa2jw <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + par[3] * acov1 + alpha * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha0 and alpha1 are fixed
# par is (log(range), alpha2)
nllh_estfsa01jw <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + alpha[2] * acov1 + par[2] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha0 and alpha2 are fixed
# par is (log(range), alpha1)
nllh_estfsa02jw <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + par[2] * acov1 + alpha[2] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha1 and alpha2 are fixed
# par is (log(range), alpha0)
nllh_estfsa12jw <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + alpha[1] * acov1 + alpha[2] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# smooth is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (log(range))
nllh_estfsa012jw <- function(par, z, lDoF, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                             pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
  par <- c(lDoF,par[1])
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# alpha0 is fixed
# par is (log(DoF), alpha1, alpha2)
nllh_estrsa0jw <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha + par[2] * acov1 + par[3] * acov2 
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# range is fixed
# smooth is fixed
# alpha1 is fixed
# par is (log(DoF), alpha0, alpha2)
nllh_estrsa1jw <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + alpha * acov1 + par[3] * acov2 
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# range is fixed
# smooth is fixed
# alpha2 is fixed
# par is (log(DoF), alpha0, alpha1)
nllh_estrsa2jw <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + par[3] * acov1 + alpha * acov2 
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# alpha0 and alpha1 are fixed
# par is (log(DoF), alpha2)
nllh_estrsa01jw <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + alpha[2] * acov1 + par[2] * acov2 
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# range is fixed
# smooth is fixed
# alpha0 and alpha2 are fixed
# par is (log(DoF), alpha1)
nllh_estrsa02jw <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + par[2] * acov1 + alpha[2] * acov2 
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# range is fixed
# smooth is fixed
# alpha1 and alpha2 are fixed
# par is (log(DoF), alpha0)
nllh_estrsa12jw <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[2] + alpha[1] * acov1 + alpha[2] * acov2 
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# range is fixed
# smooth is fixed
# alpha0, alpha1 and alpha2 are fixed
# par is (log(DoF), alpha0)
nllh_estrsa012jw <- function(par, z, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                             pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2 
  par <- c(par[1],range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

###
###

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha0 is fixed
# par is (alpha1, alpha2)
nllh_estfrsa0jw <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha + par[1] * acov1 + par[2] * acov2 
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha1 is fixed
# par is (alpha0, alpha2)
nllh_estfrsa1jw <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[1] + alpha * acov1 + par[2] * acov2 
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha2 is fixed
# par is (alpha0, alpha1)
nllh_estfrsa2jw <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                           pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par[1] + par[2] * acov1 + alpha * acov2 
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha0 and alpha1 are fixed
# par is (alpha2)
nllh_estfrsa01jw <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + alpha[2] * acov1 + par * acov2 
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed) 
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha0 and alpha2 are fixed
# par is (alpha1)
nllh_estfrsa02jw <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- alpha[1] + par * acov1 + alpha[2] * acov2 
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# lDoF = log(DoF) is fixed
# range is fixed
# smooth is fixed
# alpha1 and alpha2 are fixed
# par is (alpha0)
nllh_estfrsa12jw <- function(par, z, lDoF, range, smooth, alpha, sites, hit, acov1, acov2, jw = 2, cmat, split = FALSE, parallel = FALSE, 
                            pfun = mypmvsext, args1 = list(), args2 = list(), seed = 123) 
{
  
  Alpha <- par + alpha[1] * acov1 + alpha[2] * acov2 
  par <- c(lDoF,range)
  nllh_ESTjw(par = par, z = z, sites = sites, hit = hit, smooth = smooth, alpha = Alpha, jw = jw, cmat = cmat,
             split = split, parallel = parallel, pfun = pfun, args1 = args1, args2 = args2, seed = seed)
}

# par is (log(range), logit(smooth/2), skintercept, sk1, sk2)

# calculates V(z)
vfun_est <- function(z, sites, DoF, lambda, smooth, alpha,
                     split = FALSE, pfun = mypmvsext, args = list(), seed) 
{
  set.seed(seed)
  d <- length(z)
  qf <- function(x, A) colSums(x * (A %*% x))
  
  if(ncol(sites) != 2) stop("sites must have two columns")
  if(nrow(sites) != d) stop("incorrect number of sites")
  hmat <- as.matrix(dist(sites))
  
  # correlation model is rho(h) = exp(-(h/lambda)^smooth)
  rhomat <- exp(-(hmat/lambda)^smooth)
  rhomat2 <- sqrt(1 - rhomat*rhomat)
  
  alphastar <- numeric(d)
  for(i in 1:d) {
    rhoveci <- rhomat[i,-i]
    rhomati <- rhomat[-i,-i]
    alphastar[i] <- (alpha[i] + sum(rhoveci * alpha[-i])) /
      sqrt(1 + qf(alpha[-i], rhomati - rhoveci %o% rhoveci))
  }
  mplus <- exp((DoF/2) * log(2) - log(pi)/2 + lgamma((DoF+1)/2) + 
                 pt(alphastar * sqrt(DoF + 1), df = DoF + 1, log.p = TRUE))
  zm <- z * mplus
  
  vec <- numeric(d)
  for(i in 1:d) {
    rhoveci <- rhomat[i,-i]
    rhoveci2 <- rhomat2[i,-i]
    rhomati <- rhomat[-i,-i]
    
    eta <- sqrt(DoF+1) * ((zm[i]/zm[-i])^(-1/DoF) - rhoveci) / rhoveci2 
    ommat <- (rhomati - (rhoveci %o% rhoveci))
    # same as rmat <- ommat / (rhoveci2 %o% rhoveci2) for normalizing to correlation
    dimat <- diag(1/sqrt(diag(ommat)), nrow = nrow(ommat)) 
    rmat <- dimat %*% ommat %*% dimat
    Eigen <- tryCatch(eigen(rmat)$values, error=function(e) -1) 
    if(any(Eigen<=0)){return(-1e300)}
    
    alphaz <- sqrt(diag(ommat)) * alpha[-i] 
    extz <- sqrt(DoF + 1) * (alpha[i] + sum(rhoveci * alpha[-i]))
    vec[i] <- do.call(pfun, c(list(upper=eta, sigma=rmat, df = DoF+1, alpha = alphaz, ext = extz), args))
  }
  
  if(split) return(vec/z)
  return(sum(vec/z))
}

# calculates sum of log(-V_s(z)) over s in hitting scenario
# slst is a list containing the hitting scenario
vdfun_est <- function(z, sites, DoF, lambda, smooth, alpha, slst = list(1:d), 
                      log.p = TRUE, split = FALSE, pfun = mypmvsext, args = list(), seed) 
{
  set.seed(seed)
  qf <- function(x, A) colSums(x * (A %*% x))
  ldetinv <- function(mat) {
    cfac <- chol(mat)
    list(ldet = 2*sum(log(diag(cfac))), inv = chol2inv(cfac))
  }
  
  d <- length(z)
  if(ncol(sites) != 2) stop("sites must have two columns")
  if(nrow(sites) != d) stop("incorrect number of sites")
  tst <- sort(unlist(slst))
  if(length(tst) != d || !all(tst == 1:d)) 
    stop("invalid hitting scenario")
  
  hmat <- as.matrix(dist(sites))
  
  # correlation model is rho(h) = exp(-(h/lambda)^smooth)
  rhomat <- exp(-(hmat/lambda)^smooth)
  
  # variances are one
  sigma <- rhomat 
  isigma <- solve(sigma)
  
  alphastar <- numeric(d)
  for(i in 1:d) {
    rhoveci <- rhomat[i,-i]
    rhomati <- rhomat[-i,-i]
    alphastar[i] <- (alpha[i] + sum(rhoveci * alpha[-i])) /
      sqrt(1 + qf(alpha[-i], rhomati - rhoveci %o% rhoveci))
  }
  mplus <- exp((DoF/2) * log(2) - log(pi)/2 + lgamma((DoF+1)/2) + 
                 pt(alphastar * sqrt(DoF + 1), df = DoF + 1, log.p = TRUE))
  
  # calculates log(-V_{sv}) for log.p TRUE
  vsfun <- function(z, sites, sv = 1:d, log.p = TRUE) 
  {
    
    s <- length(sv)
    zs <- z[sv]; zo <- z[-sv]
    alphas <- alpha[sv]; alphao <- alpha[-sv] 
    mpluss <- mplus[sv]; mpluso <- mplus[-sv]
    alphastars <- alphastar[sv]; alphastaro <- alphastar[-sv]
    zms <- zs * mpluss
    zmo <- zo * mpluso
    
    sigmas <- sigma[sv, sv, drop = FALSE]
    ldi <- ldetinv(sigmas)
    ldsigmas <- ldi$ldet; isigmas <- ldi$inv
    afun <- qf(zms^(1/DoF), isigmas)
    
    sigmaso <- sigma[-sv, -sv, drop = FALSE]
    sigmaxx <- sigma[-sv, sv, drop = FALSE]
    
    if(d != s) {
      # alphabar <- as.numeric(alphas + isigmas %*% t(sigmaxx) %*% alphao /
      #                          sqrt(1 + qf(alphao, solve(sigmaso - sigmaxx %*% isigmas %*% t(sigmaxx)))))
      alphabar <- as.numeric( (alphas + isigmas %*% t(sigmaxx) %*% alphao) /
                                sqrt(1 + qf(alphao, sigmaso - sigmaxx %*% isigmas %*% t(sigmaxx))))
    } else {
      alphabar <- alpha
    }
    alphatilde <- sum(alphabar * zms^(1/DoF)) / sqrt(afun)
    
    # if sv is 1:d then mypmvt is taken as unity
    if(d != s) {
      tmp <- sigmaxx %*% isigmas 
      Sig <- sigmaso - tmp %*% t(sigmaxx)
      if(any(diag(Sig)<=0)){return(-1e300)}
      Gam <- afun / (s + DoF) * makeSymmetric(Sig)
      # Gam <- afun / (s + DoF) * makeSymmetric(sigmaso - tmp %*% t(sigmaxx))
      Eigen <- tryCatch(eigen(Gam)$values, error=function(e) -1) 
      if(any(Eigen<=0)){return(-1e300)}
      #if(!isTRUE(all.equal(Gam, t(Gam)))){   # in case of rounding issues that make the covariance matrix not a covariance matrix
      #  #warning("doing rounding of Gam matrix")
      #  Gam <- round(Gam, abs( floor( log10(max(Gam - t(Gam))) ) )-2  )
      #}
      mu <- as.numeric(tmp %*% zms^(1/DoF))
      alphavec <- alphao * sqrt((s + DoF)/afun) * sqrt(diag(Gam))
      extval <- sqrt((s + DoF)/afun) * t(alphas + t(tmp) %*% alphao) %*% zms^(1/DoF) 
      #p1 <- as.numeric(log(do.call(pfun, c(list(upper=zo^(1/DoF) - mu, sigma=Gam, df = DoF+s), args))))
      # to avoid bug in mvtnorm (reported so should be fixed soon) NOW FIXED!
      #if(length(zo)==1){
      #  p1 <- as.numeric(log(do.call(pfun, c(list(upper=as.numeric((zmo^(1/DoF) - mu)/sqrt(Gam)), sigma=as.matrix(1), df = DoF+s, alpha = alphavec, ext = extval), args))))
      #}else{
      #  p1 <- as.numeric(log(do.call(pfun, c(list(upper=zmo^(1/DoF) - mu, sigma=Gam, df = DoF+s, alpha = alphavec, ext = extval), args)))) 
      #}
      p1 <- as.numeric(log(do.call(pfun, c(list(upper=zmo^(1/DoF) - mu, sigma=Gam, df = DoF+s, alpha = alphavec, ext = extval), args))))
      
    } else {
      p1 <- 0
    }
    
    # uses DoF/2 * log(2) as 2 multiplier derives from dnorm term
    p2 <- DoF/2 * log(2) - (s-1) * log(DoF) - s/2 * log(pi) - ldsigmas/2
    p3 <- lgamma((s + DoF)/2) + (1 - DoF)/DoF * sum(log(zs)) + 1/DoF * sum(log(mpluss)) - (s + DoF)/2 * log(afun)   
    p4 <- pt(alphatilde * sqrt(s + DoF), df = s + DoF, log.p = TRUE)
    out <- p1 + p2 + p3 + p4
    if(!log.p)  out <- exp(out)
    out
  }  
  
  val <- numeric(length(slst))
  for (k in 1:length(slst)) {
    val[k] <- vsfun(z = z, sites = sites, sv = slst[[k]], log.p = TRUE) 
  }
  if(!split) val <- sum(val)
  if(!log.p) val <- exp(val)
  val
}


