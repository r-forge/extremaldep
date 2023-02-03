plot_ExtDep.np <- function(out, type, summary.mcmc, burn, y, probs, 
                           A_true, h_true, est.out, mar1, mar2, dep, 
                           QatCov1=NULL, QatCov2=QatCov1, P, 
                           labels=c(expression(y[1]),expression(y[2])), 
                           CEX=1.5, xlim, ylim, col.data, 
                           col.Qfull, col.Qfade, data=NULL, ...){
  
  # out is the output of the fExtDep.np() function
  
  est <- out$method
  ests <- c("Bayesian", "Frequentist", "Empirical")
  if(!any(est == ests)){ stop("est needs to be `Bayesian` or `Frequentist`")}
  
  if(type == "Qsets"){
    
    # Define the density function for the exponent measure:
    den <- function(w, beta, gamma){
      return(2 * dh(w, beta) * (2 * dh(w, beta) * (w)^(1-gamma[1]) * (1-w)^(1-gamma[2]) / gamma[1] / gamma[2])^(-1/(1+gamma[1]+gamma[2])))
    }
    
    # Bound function
    ghat_fun <- function(w, hhat, gamma1, gamma2){
      return((2*hhat * (w)^(1-gamma1) * (1-w)^(1-gamma2) / gamma1 / gamma2)^(1/(1+gamma1+gamma2)))
    }
  }
  
  if(est == "Bayesian"){
    
    x <- seq(from=0.0001, to=0.9999, length=100)
    nsim <- out$nsim
    
    if(missing(burn)){
      
      if(missing(summary.mcmc)){
        burn <- round(length(out$k))
        message("No burnin period provided, half of posterior discarded")
      }else{
        burn <- summary.mcmc$burn
      }
      
    }
    
    if(missing(summary.mcmc)){
      summary.mcmc <- summary_ExtDep(mcmc=out, burn=burn)
    }
    
    types <- c("summary", "returns", "A", "h", "pm", "k", "Qsets")
    if(!any(type == types)){ stop("Wrong type of graphical summary for a Bayesian estimation method")}
    
    if(type == "Qsets"){

      if(!("threshold" %in% names(out))){
        stop("Qsets only developped for raw data not maxima")
      }
      if(missing(P)){
        stop("Missing probabilities (argument `P')")
      }
      if(length(P)>3){
        stop("Length of `P' cannot be greater than 3")
      }
      kn <- out$kn

      if(dep==TRUE && missing(data)){
        stop("data must be provided when dependence is plotted")
      }
      
    }
    
    ###################### START PLOT RETURNS ####################
    
    ################################################################################
    # INPUT:                                                                     ###
    # y is a (m x 2)-dimensional matrix of the thresholds                        ###
    # probs is the output obtained with returns function, i.e. vector of       ###
    #          returns values                                                    ###
    ################################################################################ 
    
    plot.returns <- function(y, probs, CEX=1.5,
                             labels=c(expression(y[1]),expression(y[2])), data=NULL, ...){
      
      op1 <- par(mar = c(1, 1, 0, 0), oma = c(3, 4, 0.5, 0.5), mgp = c(1, 1, 0), cex.axis=CEX)
      on.exit(par(op1))
      
      if(is.vector(y)){
        plot(y,probs,type='l',col=2,lwd=2.5,ylim=range(probs,na.rm=T), ylab="", xlab="", ...)
        mtext('y',side=1,line=3,cex=CEX)
        mtext(expression(P(Y[1]>y,Y[2]>y)),side=2,line=3,cex=CEX)  
      }
      else{
        ny <- sqrt(dim(y)[1])
        yy <- y[1:ny,1]
        col <- gray(100:50/100)
        plot(yy, yy, type='n', ylab='', xlab='', ...)
        image(yy, yy, matrix(probs, ny), xlab='', ylab='', col=col, add=T)
        contour(yy, yy, matrix(probs,ny), add=T,col=2,lwd=2, labcex=CEX)
        if(!is.null(data)){
          points(data, pch=16)
        }
        mtext(labels[1],side=1,line=3,cex=CEX)
        mtext(labels[2],side=2,line=3,cex=CEX)
      }
    }
    ###################### END PLOT RETURNS ####################
    
    ###################### START PLOT PICKANDS ####################
    
    ################################################################################
    # INPUT:                                                                     ###
    # w is a bidimensional unit-simplex                                          ###
    # summary.mcmc is the output obtained with summary.bbeed function            ###
    ################################################################################ 
    
    plot.A <- function(w, summary.mcmc, CEX=1.5, ...){
      # summary.mcmc = output PostMCMC
      op3 <- par(cex.axis=CEX)
      on.exit(par(op3))
      
      A_hat <- summary.mcmc$A.mean
      lowA <- summary.mcmc$A.low
      upA <- summary.mcmc$A.up
      
      plot(w, A_hat, type="n", xlim=c(0,1), 
           ylim=c(.5, 1), ylab="", xlab="", ...)
      polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = 1, lwd = 1, border = 'grey')
      polygon(c(rev(w), w), c(rev(lowA), upA), col = 'grey80', border = NA)
      lines(w, A_hat, lwd=2, col=2)
      
      mtext('t',side=1,line=3,cex=CEX)
      mtext('A(t)',side=2,line=3,cex=CEX)
    }
    ###################### END PLOT PICKANDS ####################
    
    # Median curve Pickands + Credibility bands
    # library(fda)
    
    fbplot_A <- function(A_post, A_true, wgrid, CEX=1.5){
      # A_post = PostMCMC$A_post
      # is the matrix of Pickands from the final chains of the MCMC
      # which refer to those k in (q1, q3) of the posterior
      # A_true = True Pickands
      op4 <- par(oma = c(3, 4, 0.5, 0.5),mgp = c(1, 1, 0),cex.axis=CEX)
      on.exit(par(op4))
      
      A_med <- fbplot(t(A_post),wgrid,xlim=c(0,1),ylim=c(0.5,1),
                      method='BD2',color='grey70',xlab='',ylab='',
                      outliercol='grey50',barcol='grey80') 
      lines(wgrid,A_true,col=1,lwd=2)
      lines(wgrid,A_post[A_med$medcurve[1],],col=2,lwd=2)
      polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = 1, lwd = 1, border = 'grey')
      mtext('w',side=1,line=3,cex=CEX)
      mtext('A(w)',side=2,line=3,cex=CEX)
    }
    
    
    ###################### START PLOT ANGULAR DISTRIBUTION ####################
    
    ################################################################################
    # INPUT:                                                                     ###
    # w is a bidimensional unit-simplex                                          ###
    # summary.mcmc is the output obtained with summary.bbeed function            ###
    ################################################################################ 
    
    plot.h <- function(w, summary.mcmc, CEX=1.5, ...){ 
      
      # summary.mcmc = output summary.mcmc
      # pm = theorethical point masses, zero by default
      op5 <- par(cex.axis=CEX)
      on.exit(par(op5))
      
      nw <- length(w)
      
      h_hat <- summary.mcmc$h.mean
      lowh <- summary.mcmc$h.low
      uph <- summary.mcmc$h.up
      
      p0_hat <- summary.mcmc$p0.mean
      lowp0 <- summary.mcmc$p0.low
      upp0 <- summary.mcmc$p0.up
      
      p1_hat <- summary.mcmc$p1.mean
      lowp1 <- summary.mcmc$p1.low
      upp1 <- summary.mcmc$p1.up
      
      ylim=c(0, max(uph, h_hat, na.rm=T))
      
      plot(w, h_hat, type="n", xlim=c(0,1), ylim=ylim, ylab="", xlab="",...)
      
      polygon(c(rev(w), w), c(rev(lowh), uph), col = 'grey80', border = NA)
      points(0, lowp0 , pch=16, cex=2,col='grey80')
      points(0, upp0 , pch=16, cex=2,col='grey80')
      points(1, lowp1 , pch=16, cex=2,col='grey80')
      points(1, upp1 , pch=16, cex=2,col='grey80')
      
      lines(w[-c(1,nw)], h_hat[-c(1,nw)], lwd=2, col=2)
      points(0, p0_hat , pch=16, cex=2,col=2)
      points(1, p1_hat , pch=16, cex=2,col=2)
      
      mtext(expression(omega),side=1,line=3,cex=CEX)
      mtext(expression(h(omega)),side=2,line=3,cex=CEX)  
    }
    
    ###################### END PLOT ANGULAR DISTRIBUTION ####################
    
    fbplot_h <- function(pmcmc, h_true, wgrid, CEX=1.5){
      # A_post = PostMCMC$A_post
      # is the matrix of Pickands from the final chains of the MCMC
      # which refer to those k in (q1, q3) of the posterior
      # A_true = True Pickands
      op6 <- par(oma = c(3, 4, 0.5, 0.5), mgp = c(1, 1, 0),cex.axis=CEX)
      on.exit(par(op6))
      
      h_med <- fbplot(t(pmcmc$h_post),wgrid,xlim=c(0,1),ylim=c(0,max(h_true)+.5),
                      method='BD2',color='grey80',xlab='',ylab='',
                      outliercol='white',barcol='white',fullout=F) 
      lines(wgrid,pmcmc$h_post[h_med$medcurve[1],],col='grey80',lwd=4)
      lines(wgrid,h_true,col=1,lwd=2)
      lines(wgrid,pmcmc$h.mean,col=2,lwd=2)
      mtext('u',side=1,line=3,cex=CEX)
      mtext('h(u)',side=2,line=3,cex=CEX)
    }
    
    ###################### START PLOT PRIOR VS POSTERIOR k ####################
    
    ################################################################################
    # INPUT:                                                                     ###
    # mcmc is the output obtained with bbeed function                            ###
    # burn is the number of samples to discard                                   ###
    # nsim is the number of the iterations of the chain                          ###
    ################################################################################ 
    
    PriorVSPosterior.k <- function(mcmc, nsim, burn, CEX=1.5, ...){
      # mcmc = output bbeed
      
      op7 <- par(cex.axis=CEX)
      on.exit(par(op7))
      
      prior.k <- mcmc$prior$k
      
      chain.k <- mcmc$k[burn:nsim]
      t <- table(chain.k)
      kval <- as.numeric(names(t))
      
      if(prior.k=='pois') priork <- dpois(0:max(kval)-3,mcmc$prior$hyperparam$mu)
      if(prior.k=='nbinom') priork <- dnbinom(0:max(kval)-3,size=mcmc$prior$hyperparam$mu.nbinom^2 / (mcmc$prior$hyperparam$var.nbinom-mcmc$prior$hyperparam$mu.nbinom),prob=mcmc$prior$hyperparam$mu.nbinom/mcmc$prior$hyperparam$var.nbinom)
      
      xlim <- c(min(chain.k), max(chain.k))
      postk <- as.vector(t)/length(chain.k)
      ylim <- c(0,max(postk,priork))
      
      plot(0:max(kval),priork,type="h",xlim=xlim,ylim=ylim,
           lwd=3,col="chartreuse3",ylab="",xlab="",...)
      points(0:max(kval),priork,pch=16,cex=2,col="green2")
      for(i in 1:length(kval))
        lines( c(kval[i],kval[i]), c(0,postk[i]), col = "dark red", lwd = 2)
      points(kval,postk,pch=16,cex=2,col="red")
      mtext(expression(kappa),side=1,line=3,cex=CEX)
    }
    
    ###################### END PLOT PRIOR VS POSTERIOR k ####################
    
    ###################### START PLOT PRIOR VS POSTERIOR p0 ####################
    
    ################################################################################
    # INPUT:                                                                     ###
    # mcmc is the output obtained with bbeed function                            ###
    # burn is the number of samples to discard                                   ###
    # nsim is the number of the iterations of the chain                          ###
    ################################################################################ 
    
    PriorVSPosterior.pm <- function(mcmc, nsim, burn, CEX=1.5, ...){
      
      op8 <-par(cex.axis=CEX)
      on.exit(par(op8))
      
      prior.pm <- mcmc$prior$pm
      chain.p0 <- mcmc$pm[burn:nsim,1]
      
      postp0 <- hist(chain.p0, plot=F)$density
      ydmax <- max(postp0,2)
      
      a <- mcmc$prior$hyperparam$a
      b <- mcmc$prior$hyperparam$b
      
      if(prior.pm == 'unif'){
        if(length(a)>1) stop('Check hyperparameters a and b.')
        curve(dunif(x,a,b),0,.5,ylim=c(0,ydmax+1),col='green2',lwd=2, xlab="", ylab="", ...)
      }
      
      if(prior.pm == 'beta'){
        if(length(a)==1) stop('Check hyperparameters a and b.')
        curve(dbeta(x,a,b),0,.5,ylim=c(0,ydmax+1),col='green2',lwd=2, xlab="", ylab="", ...)
      }
      
      lines(seq(0,.5,length=length(postp0)), postp0, col = "red", lwd = 2)
      mtext(expression(p[0]),side=1,line=3,cex=CEX)
    }
    
    ###################### START PLOT EXTREME QUANTILE SETS ########################
    
    ################################################################################
    # INPUT:                                                                     ###
    # mcmc is the output obtained with fExtDep.np() function                     ###
    # burn is the number of samples to discard                                   ###
    # nsim is the number of the iterations of the chain                          ###
    ################################################################################ 
    
    Plot.Qset <- function(mcmc, nsim, burn, kn,  summary.mcmc, est.out,
                          QatCov1=NULL, QatCov2=NULL, P, 
                          mar1, mar2, dep=TRUE,  
                          xlim, ylim, xlab=NULL, ylab=NULL, CEX=1.5, 
                          col.data, col.Qfull, col.Qfade,...){
      
      # preliminary function
      
      func <- function(y){
        y <- y[is.finite(y)]
        return(c(quantile(y, 0.05), mean(y), quantile(y, 0.95)))
      }
      
      mar.fit <- mcmc$mar.fit
      
      nw <- 100
      w <- seq(0.00001, .99999, length=nw)
      m <- length(P)
      Qhat <- Qhat_up <- Qhat_low <- array(0, c(nw,2,m))
      npost <- nsim-burn+1
      
      if(mar.fit){
        d.cov1 <- ncol(mcmc$mar1) -2 
        d.cov2 <- ncol(mcmc$mar2) -2 
        
        if(is.null(QatCov1) && d.cov1>1){stop("QatCov1 argument required")}
        if(is.null(QatCov2) && d.cov2>1){stop("QatCov2 argument required")}
        if(!is.null(QatCov1) && (d.cov1>1) && (ncol(QatCov1) != (d.cov1-1))){stop("Wrong dimensions of QatCov1")}
        if(!is.null(QatCov2) && (d.cov2>1) && (ncol(QatCov2) != (d.cov2-1))){stop("Wrong dimensions of QatCov2")}
        
        if(d.cov1==1){
          n.cov1 <- 1
          Cov1 <- as.matrix(1)
        }else{
          n.cov1 <- nrow(QatCov1)
          Cov1 <- cbind(rep(1, n.cov1),QatCov1)
        }
        if(d.cov2==1){
          n.cov2 <- 1
          Cov2 <- as.matrix(1)
        }else{
          n.cov2 <- nrow(QatCov2)
          Cov2 <- cbind(rep(1, n.cov2),QatCov2)
        }
        
        if(n.cov1 != n.cov2){stop("QatCov1 and QatCov2 sould have the same dimensions")}
        if(n.cov1>3){stop("plot.ExtDep will display maximum 3 covariate levels")}
        
        muhat1_post <- summary.mcmc$mar1_post[,1:(ncol(Cov1))] %*% t(Cov1) # A npost by nrow(Cov1) matrix
        muhat2_post <- summary.mcmc$mar2_post[,1:(ncol(Cov2))] %*% t(Cov2) # A npost by nrow(Cov2) matrix
        sighat1_post <- summary.mcmc$mar1_post[,ncol(Cov1)+1]
        sighat2_post <- summary.mcmc$mar2_post[,ncol(Cov2)+1]  
        gamhat1_post <- summary.mcmc$mar1_post[,ncol(Cov1)+2]
        gamhat2_post <- summary.mcmc$mar2_post[,ncol(Cov2)+2]
        
      }else{
        if(missing(mar1)){stop("mar1 should be provided when marginal parameters weren't fitted")}
        if(missing(mar2)){stop("mar2 should be provided when marginal parameters weren't fitted")}
        muhat1_post <- matrix(rep(mar1[1], npost), ncol=1)
        muhat2_post <- matrix(rep(mar2[1], npost), ncol=1)
        sighat1_post <- rep(mar1[2], npost)
        sighat2_post <- rep(mar2[2], npost)
        gamhat1_post <- rep(mar1[3], npost)
        gamhat2_post <- rep(mar2[3], npost)
        Cov1 <- Cov2 <- as.matrix(1)
        n.cov1 <- n.cov2 <- 1
      }

      hhat <- rbind(summary.mcmc$h.low,summary.mcmc$h.mean,summary.mcmc$h.up)
      
      # # Define the density function for the exponent measure:
      # den <- function(w, beta, gamma){
      #   return(2 * dh(w, beta) * (2 * dh(w, beta) * (w)^(1-gamma[1]) * (1-w)^(1-gamma[2]) / gamma[1] / gamma[2])^(-1/(1+gamma[1]+gamma[2])))
      # }
      # 
      # # Bound function
      # ghat_fun <- function(w, hhat, gamma1, gamma2){
      #   return((2*hhat * (w)^(1-gamma1) * (1-w)^(1-gamma2) / gamma1 / gamma2)^(1/(1+gamma1+gamma2)))
      # }
      
      ###############################################################
      # START Estimation of the basic-set S and measure nu(S) :
      ###############################################################
      
      if(missing(est.out)){
        
        # Estimation of the bound
        ghat_post <- matrix(nrow = npost, ncol = nw)
        for(i in 1:nw) {
          ghat_post[,i] <- ghat_fun(w[i], hhat=summary.mcmc$h_post[,i], gamma1 = gamhat1_post, gamma2 = gamhat2_post)
        }                                                                                                                        
        ghat <- apply(ghat_post,2, func)
        
        # Estimation of the basic-set S:
        Shat_post <- array(0, c(nw,2,npost))
        for(j in 1:npost) Shat_post[,,j] <- cbind(ghat_post[j,]*w, ghat_post[j,]*(1-w)) 
        
        Shat <- array(0, c(nw,2,3))
        for(j in 1:3) Shat[,,j] <- cbind(ghat[j,]*w, ghat[j,]*(1-w)) 
        
        # Estimation of the measure nu(S):
        nuShat_post <- numeric(npost)
        nuShat <- NULL
        
        for(i in 1:npost) nuShat_post[i] <- integrate(Vectorize(function(t){den(t,beta=summary.mcmc$beta_post[[i]], gamma=c(gamhat1_post[i],gamhat2_post[i]))}), 0,1)$value
        nuShat <- func(nuShat_post)
        
        est.out <- list(ghat=ghat, Shat=Shat, Shat_post=Shat_post, nuShat=nuShat, nuShat_post=nuShat_post)  
      }else{
        ghat <- est.out$ghat 
        Shat <- est.out$Shat 
        Shat_post <- est.out$Shat_post
        nuShat <- est.out$nuShat
        nuShat_post <- est.out$nuShat_post
      }
      
      # Graphical representations
      
      if(dep){
        par(mfrow=c(2,3), mar=c(3, 3, 0.5, 0.2), mgp=c(2,.8,0))
      }else{
        par(mfrow=c(1,m), mar=c(3, 3, 0.5, 0.2), mgp=c(2,.8,0))
      }
      
      if(is.null(col.data)){ 
        rbPal <- colorRampPalette(c('blue','red'))
      }  
      
      # Illustrations - Dependence
      if(dep){
        
        plot(NULL, xlab="w", ylab=expression(1/q[symbol("\052")](w)), ylim=c(0,max(ghat)), xlim=c(0,1))
        polygon(c(w, rev(w)), c(ghat[3,], rev(ghat[1,])), col="gray")
        lines(w, ghat[2,], lwd=2, lty=3)
        
        plot(NULL, xlim=c(0,max(Shat[,1,])), ylim=c(0,max(Shat[,2,])), xlab=xlab, ylab=ylab)
        polygon(c(Shat[,1,3], rev(Shat[,1,1])), c(Shat[,2,3], rev(Shat[,2,1])), col="gray")
        points(Shat[,,2], type="l", col="gray0", lwd=2, lty=3)
        
        if(is.null(col.data)){ 
          if(n.cov1 > 1){
            Colors <- rbPal(10)
            col.data <- Colors[as.numeric(cut(out$cov1[,2],breaks = 10))]
          }else{
            col.data <- "black" #rbPal(1)
          }
        }  
        
        plot(data, pch=19, xlab=xlab, ylab=ylab, xlim=c(0,max(data[,1])), ylim=c(0,max(data[,2])), col=col.data )
        
      }
      
      ###############################################################
      # START Estimation of the Extreme quantile sets :
      ###############################################################
      
      if(is.null(col.Qfull) || is.null(col.Qfade)){
        col.Qfull <- adjustcolor( rbPal(n.cov1), alpha.f = 1)
        col.Qfade <- adjustcolor( rbPal(n.cov1), alpha.f = 0.6)
      }  
      
      if(is.null(xlim)) xlim <- c(0,max(data[,1]))
      if(is.null(ylim)) ylim <- c(0,max(data[,2]))
      
      # Computation of the quantiles
      for(i in 1:m){
        
        plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
        
        for(z in 1:n.cov1){
          
          tmp <- array(0, c(nw,2,npost))
          tmp2 <- array(0, c(nw,2,3)) # 3 because we compute the mean and 90% cred. regions
          
          for(j in 1:nw){
            
            tmp <- rbind( muhat1_post[,z] + sighat1_post / gamhat1_post  * (( kn[1] * nuShat_post*Shat_post[j,1,]/P[i])^(gamhat1_post) -1) ,
                          muhat2_post[,z] + sighat2_post / gamhat2_post  * (( kn[2] * nuShat_post*Shat_post[j,2,]/P[i])^(gamhat2_post) -1) )
            
            tmp2[j,,] <- cbind(c(tmp[1,which(tmp[1,]==func(tmp[1,])[1])[1]],tmp[2,which(tmp[2,]==func(tmp[2,])[1])[1]]), 
                               apply(tmp,1,mean),
                               c(tmp[1,which(tmp[1,]==func(tmp[1,])[3])[1]],tmp[2,which(tmp[2,]==func(tmp[2,])[3])[1]])) 
          }
          
          polygon( c(tmp2[-1,1,3], rev(tmp2[-1,1,1])),  c(tmp2[-1,2,3], rev(tmp2[-1,2,1])), col=col.Qfade[z], border=NA )
          lines(tmp2[-1,,2], col=col.Qfull[z], lwd=2)
          
          assign(paste("Qset_P",i,"_CovNum_",z,"_post",sep=""), tmp)
          assign(paste("Qset_P",i,"_CovNum_",z,sep=""), tmp2) # 3 because we compute the mean and 90% cred. regions
          
        }
        
      }
      
      q.out <- mget(ls(pattern = "Qset_P"))
      
      return(list(est.out=est.out, q.out=q.out))
    }  
    
    ###################### END PLOT EXTREME QUANTILE SETS ##########################  
      
    if(type == "returns"){
      probs <- returns(out=out, summary.mcmc=summary.mcmc, y=y)
      plot.returns(y=y, probs=probs, CEX=CEX, labels=labels, data=data, ...)
    } 
    if(type == "A"){
      if(!missing(A_true)){
        fbplot_A(A_post=summary.mcmc$A_post, A_true, wgrid=summary.mcmc$w, CEX=CEX)
      }else{
        plot.A(w=x, summary.mcmc=summary.mcmc, CEX=CEX, ...)
      }
    }
    if(type == "h"){
      if(!missing(h_true)){
        fbplot_h(pmcmc=summary.mcmc, h_true, wgrid=summary.mcmc$w, CEX=CEX)
      }else{
        plot.h(w=x, summary.mcmc=summary.mcmc, CEX=CEX, ...)
      }
    }
    if(type == "pm"){
      PriorVSPosterior.pm(mcmc=out, nsim=nsim, burn=burn, CEX=CEX, ...)
    }
    if(type == "k"){
      PriorVSPosterior.k(mcmc=out, nsim=nsim, burn=burn, CEX=CEX, ...)
    }
    if(type == "summary"){
      oldpar1 <- par(mfrow=c(2,2), pty='s', mar = c(4.5, 4.2, 0.25, 0.25), 
                     oma = c(0, 0, 0, 0), mgp = c(1, 1, 0), cex.axis=CEX)
      on.exit(par(oldpar1))
      
      plot.A(w=x, summary.mcmc=summary.mcmc, ...)
      PriorVSPosterior.k(mcmc=out, nsim=nsim, burn=burn, ...)
      plot.h(w=x, summary.mcmc=summary.mcmc, ...)
      PriorVSPosterior.pm(mcmc=out, nsim=nsim, burn=burn, ...)
    }  
    
    if(type == "Qsets"){
      
      if(missing(col.data)) col.data <- NULL
      if(missing(col.Qfull)) col.Qfull <- NULL
      if(missing(col.Qfade)) col.Qfade <- NULL
      if(missing(xlim)) xlim <- NULL
      if(missing(ylim)) ylim <- NULL

      Plot.Qset(mcmc=out, nsim=nsim, burn=burn, kn=kn, summary.mcmc=summary.mcmc,
                est.out=est.out, QatCov1=QatCov1, QatCov2=QatCov2,
                P=P, mar1=mar1, mar2=mar2, dep=dep,  
                xlim=xlim, ylim=ylim, xlab=labels[1], ylab=labels[2], CEX=1.5, 
                col.data=col.data, col.Qfull=col.Qfull, col.Qfade=col.Qfade,...)
      
    }
    
  }else if(est == "Frequentist"){
    
    if(type %in% c("returns", "pm", "k")){stop("Wrong type specified when est='Frequentist'")}
    
    if(type == "Qsets"){
      
      if(out$type == "maxima"){
        stop("Qsets only developped for raw data not maxima")
      }
      if(missing(P)){
        stop("Missing probabilities (argument `P')")
      }
      if(length(P)>3){
        stop("Length of `P' cannot be greater than 3")
      }
      kn <- out$kn
      
      if(dep==TRUE && missing(data)){
        stop("data must be provided when dependence is plotted")
      }
      
    }
    
    plot.A.F <- function(out, CEX=1.5, ...){
  
      op3 <- par(cex.axis=CEX)
      on.exit(par(op3))
      
      plot(out$w, out$Ahat$A, type="l", xlim=c(0,1), 
           ylim=c(.5, 1), ylab="", xlab="", lwd=2, col=2, ...)
      polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = 1, lwd = 1, border = 'grey')

      mtext('t',side=1,line=3,cex=CEX)
      mtext('A(t)',side=2,line=3,cex=CEX)
    }
    
    plot.h.F <- function(out, CEX=1.5, ...){ 
      
      # summary.mcmc = output summary.mcmc
      # pm = theorethical point masses, zero by default
      op5 <- par(cex.axis=CEX)
      on.exit(par(op5))

      ylim=c(0, max(out$hhat, na.rm=T))
      
      hist(out$extdata[,2]/rowSums(out$extdata), prob=TRUE, 
           xlim=c(0,1), main="", ylim=ylim, ylab="", xlab="",...)
      
      lines(out$w, out$hhat, type="l", lwd=2, col=2)
      
      points(0, head(out$p0,1), pch=16, cex=2,col=2)
      points(1, tail(out$p1,1), pch=16, cex=2,col=2)
      
      mtext(expression(omega),side=1,line=3,cex=CEX)
      mtext(expression(h(omega)),side=2,line=3,cex=CEX)  
    }
    
    Plot.Qset.F <- function(obj, mar1, mar2, dep=TRUE, est.out=FALSE, P,   
                          xlim=xlim, ylim=ylim, xlab=NULL, ylab=NULL, 
                          col.data, col.Qfull, CEX=1.5, ...){
      
      if(all(c("f1","f2") %in% names(obj))){ 
        mar.fit <- TRUE 
      }else{
        mar.fit <- FALSE
      }
      if(mar.fit==FALSE && (is.null(mar1) || is.null(mar2)) ){stop("missing marginal parameters mar1 and mar2")}
      
      m <- length(P)
      
      hhat <- obj$hhat
      Ahat <- obj$Ahat
      w <- obj$w
      nw <- length(w)
      kn <- 1-obj$q
      Qhat <- array(0, c(nw,2,m))
      
      if(mar.fit){
        muhat1 <- obj$f1$est[1]
        sighat1 <- obj$f1$est[2]
        gamhat1 <- obj$f1$est[3]
        muhat2 <- obj$f2$est[1]
        sighat2 <- obj$f2$est[2]
        gamhat2 <- obj$f2$est[3]  
      }else{
        muhat1 <- mar1[1]
        sighat1 <- mar1[2]
        gamhat1 <- mar1[3]
        muhat2 <- mar2[1]
        sighat2 <- mar2[2]
        gamhat2 <- mar2[3]  
      }
      
      # START Estimation of the basic-set S and measure nu(S) :
      
      if(missing(est.out)){
      
        # Estimation of the bound
        ghat <- ghat_fun(w=w, hhat=hhat, gamma1=gamhat1, gamma2=gamhat2)
        
        # Estimation of the basic-set S:
        xS <- ghat*w
        yS <- ghat*(1-w)
        Shat <- cbind(xS,yS)
        
        #Estimation of the measure nu(S):
        nuShat <- adaptIntegrate(den, 0 , 1, beta=Ahat$beta, gamma=c(gamhat1,gamhat2))$integral  
        
      }else{
        ghat <- est.out$ghat
        Shat <- est.out$Shat
        nuShat <- est.out$nuShat
      }
      
      est.out <- list(ghat=ghat, Shat=Shat, nuShat=nuShat)  
      
      # Graphical representations
      
      if(dep){
        par(mfrow=c(1,3), mar=c(3, 3, 0.5, 0.2), mgp=c(2,.8,0))
      }else{
        par(mfrow=c(1,1), mar=c(3, 3, 0.5, 0.2), mgp=c(2,.8,0))
      }
      
      if(is.null(col.data)){ 
        col.data <- "black"
      }  
      
      # Illustrations - Dependence
      if(dep){
        
        plot(w, ghat, type="l", lwd=2, xlab="w", ylab=expression(1/q[symbol("\052")](w)), ylim=c(0,max(ghat)), xlim=c(0,1))
        
        plot(Shat, type="l", lwd=2, xlim=c(0,max(Shat[,1])), ylim=c(0,max(Shat[,2])), xlab=xlab, ylab=ylab)
        
      }
      
      # START Estimation of quantile regions:
      
      if(is.null(col.Qfull)){
        col.Qfull <- rep("black", m)
      }  
      
      if(is.null(xlim)) xlim <- c(0,max(data[,1]))
      if(is.null(ylim)) ylim <- c(0,max(data[,2]))
      
      plot(data, pch=19, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,, col=col.data, ... )
      
      for(i in 1:m){
        # Qhat[,,i] <- cbind(c1*(((nuShat*xS/P[i]))^gamhat1),
        #                    c2*(((nuShat*yS/P[i]))^gamhat2))  
        Qhat[,,i] <- cbind(muhat1 + sighat1 / gamhat1 * ((kn*nuShat*xS/P[i])^gamhat1-1),
                           muhat2 + sighat2 / gamhat2 * ((kn*nuShat*yS/P[i])^gamhat2-1))  
        
        lines(Qhat[,,i], col=col.Qfull[i], lwd=2)
      }
      
      return(list(est.out=est.out, q.out=Qhat))
      
    }
    
    if(type == "A"){
      plot.A.F(out=out, CEX=CEX, ...)
    }
    if(type == "h"){
      plot.h.F(out=out, CEX=CEX, ...)
    }
    if(type == "summary"){
      oldpar1 <- par(mfrow=c(1,2), pty='s', mar = c(2, 4.25, 0.25, 0.25), 
                     oma = c(0, 0, 0, 0), mgp = c(1, 1, 0), cex.axis=CEX)
      on.exit(par(oldpar1))
      
      plot.A.F(out=out, CEX=CEX, ...)
      plot.h.F(out=out, CEX=CEX, ...)
    } 
    if(type == "Qsets"){
      
      if(missing(col.data)) col.data <- NULL
      if(missing(col.Qfull)) col.Qfull <- NULL
      if(missing(xlim)) xlim <- NULL
      if(missing(ylim)) ylim <- NULL
      
      Plot.Qset.F(obj=out, mar1=mar1, mar2=mar2, dep=dep, est.out=est.out, P=P,   
                  xlim=xlim, ylim=ylim, xlab=labels[1], ylab=labels[2], 
                  col.data=col.data, col.Qfull=col.Qfull, CEX=1.5, ...)
    }
    
    
  }else if(est == "Empirical"){
    
    if(type %in% c("returns", "pm", "k")){stop("Wrong type specified when est='Empirical'")}
    
    plot.psi <- function(out, CEX=1.5, ...){
      
      op3 <- par(cex.axis=CEX)
      on.exit(par(op3))
      
      plot(out$psi_hat, type="l", xlim=c(0,pi/2), 
           ylim=c(0, max(out$psi_hat[,2])), ylab="", xlab="", lwd=2, col=2,, xaxt = "n", ...)
      axis(side = 1, at = c((0:4)*pi/8), labels = c("0", "pi/8", "pi/4", "3*pi/8", "pi/2"))

      mtext(expression(theta),side=1,line=3,cex=CEX)
      mtext(expression(psi(theta)),side=2,line=3,cex=CEX)
    }
    
    plot.h.E <- function(out, CEX=1.5, ...){ 
      
      op5 <- par(cex.axis=CEX)
      on.exit(par(op5))
      
      ylim=c(0, max(out$h_hat[,2], na.rm=T))
      
      plot(out$h_hat, type="l", xlim=c(0,1), ylim=ylim, 
           ylab="", xlab="", lwd=2, col=2, ...)
      
      mtext(expression(omega),side=1,line=3,cex=CEX)
      mtext(expression(h(omega)),side=2,line=3,cex=CEX)  
    }
    
    Plot.Qset.E <- function(obj, mar1, mar2, dep=TRUE, est.out=FALSE, P,   
                          xlim=xlim, ylim=ylim, xlab=NULL, ylab=NULL, 
                          col.data, col.Qfull, CEX=1.5, ...){
      
      muhat1 <- mar1[1]
      sighat1 <- mar1[2]
      gamhat1 <- mar1[3]
      muhat2 <- mar2[1]
      sighat2 <- mar2[2]
      gamhat2 <- mar2[3] 
      
      m <- length(P)
      
      if(missing(est.out)){
        fi <- obj$fi
        
        w <- fi$weights
        loc <- fi$angles
        lok <- pi/2 - loc
        M <- length(loc)
        h <- 0.45
        
        ## Kernels
        kernK = function(x){if(x >= -1 && x <= 1){ret = (15/16)*((1-x^2)^2)}
          if(x < -1 || x > 1){ret = 0}
          ret}
        
        kernK1 = function(x){if(x >= -1 && x <= 1){ret = x*(15/16)*((1-x^2)^2)}
          if(x < -1 || x > 1){ret = 0}
          ret}
        
        kernK2 = function(x){if(x >= -1 && x <= 1){ret = (x^2)*(15/16)*((1-x^2)^2)}
          if(x < -1 || x > 1){ret = 0}
          ret}
        
        kernK = Vectorize(kernK)
        kernK1 = Vectorize(kernK1)
        kernK2 = Vectorize(kernK2)
        
        a0 = function(y){(15/16)*(y^5/5 - 2*y^3/3 + y + 8/15)}
        a1 = function(y){(15/16)*(y^6/6 -   y^4/2 + y^2/2 - 1/6)}
        a2 = function(y){(15/16)*(y^7/7 - 2*y^5/5 + y^3/3 + 8/105)}
        a0 = Vectorize(a0)
        a1 = Vectorize(a1)
        a2 = Vectorize(a2)
        
        kernBL = function(y){
          ly = length(y)    
          ret = c(1:ly)
          p = (y+loc/h)[1]
          i = 1
          while(i <= ly){       
            ret[i] = (  (  a2(p)-a1(p)*y[i]   )*kernK(y[i])  )/( a0(p)*a2(p)-(a1(p))^2 )
            i = i+1  } # end of while 
          ret} # end of function
        
        kernBR = function(y){
          ly = length(y)
          ret = c(1:ly)
          p = ((y*h+lok)/h)[1]
          i = 1
          while(i<=ly){
            ret[i] = (  (  a2(p)-(a1(p))*y[i]  )*kernK(y[i])  )/( a0(p)*a2(p)-(a1(p))^2 )
            i = i+1  } # end of while 
          ret} # end of function
        
        fi_hat = function(x){
          if(x >= h && x <= pi/2 - h){ret = sum((w/h)*kernK((x-loc)/h))}
          if(x >= 0 && x < h){ret = sum((w/h)*kernBL((x-loc)/h))}
          if(x >= pi/2 - h && x <= pi/2){ret = sum((w/h)*kernBR((pi/2-x-lok)/h))}
          ret}
        fi_hat = Vectorize(fi_hat)
        
        g = function(theta){
          (  ((fi_hat(theta)  )*((cos(min(theta,pi/2-h)))^(1-gamhat1)  )*(  (sin(max(theta,h)))^(1-gamhat2)  ))/(gamhat1*gamhat2))^(-1/(gamhat1+gamhat2+1))
        }
        
        ## \nu(S)
        nuS_integrand <- function(x){g(x)*fi_hat(x)}
        nuS_integrand <- Vectorize(nuS_integrand)
        nuShat <- integrate(nuS_integrand, 0, pi/2)$value
        
        ## S
        S_r <- function(theta){ 1/g(theta)}
        S_r <- Vectorize(S_r)
        
        theta <- obj$theta.seq
        Srr <- S_r(theta)
        
        xS <- Srr*cos(theta)
        yS <- Srr*sin(theta)
        Shat <- cbind(xS,yS)
        
        est.out <- list(Srr=Srr, Shat=Shat, nuShat=nuShat)  
        
      }else{
        theta <- obj$theta.seq
        Srr <- est.out$Srr
        Shat <- est.out$Shat
        nuShat <- est.out$nuShat
        xS <- Shat[,1]
        yS <- Shat[,2]
      }
      
      kn1 <- 0.1
      kn2 <- 0.1
      
      # Graphical representations
      
      if(dep){
        par(mfrow=c(1,3), mar=c(3, 3, 0.5, 0.2), mgp=c(2,.8,0))
      }else{
        par(mfrow=c(1,1), mar=c(3, 3, 0.5, 0.2), mgp=c(2,.8,0))
      }
      
      if(is.null(col.data)){ 
        col.data <- "black"
      }  
      
      # Illustrations - Dependence
      if(dep){
        
        plot(theta, Srr, type="l", lwd=2, xlab=expression(theta), ylab=expression(1/q[symbol("\052")](theta)), ylim=c(0,max(Srr)), xlim=c(0,1))
        
        plot(Shat, type="l", lwd=2, xlim=c(0,max(Shat[,1])), ylim=c(0,max(Shat[,2])), xlab=xlab, ylab=ylab)
        
      }
      
      # START Estimation of quantile regions:
      
      Qhat <- array(0, c(nrow(Shat),2,m))
      
      if(is.null(col.Qfull)){
        col.Qfull <- rep("black", m)
      }  
      
      if(is.null(xlim)) xlim <- c(0,max(data[,1]))
      if(is.null(ylim)) ylim <- c(0,max(data[,2]))
      
      plot(data, pch=19, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,, col=col.data, ... )
      
      for(i in 1:m){

        Qhat[,,i] <- cbind(muhat1 + sighat1 / gamhat1 * ((kn1*nuShat*xS/P[i])^gamhat1-1),
                           muhat2 + sighat2 / gamhat2 * ((kn2*nuShat*yS/P[i])^gamhat2-1))
        
        lines(Qhat[,,i], col=col.Qfull[i], lwd=2)
      }
      
      return(list(est.out=est.out, q.out=Qhat))
      
    }  
    
    if(type == "psi"){
      plot.psi(out=out, CEX=CEX, ...)
    }
    if(type == "h"){
      plot.h.E(out=out, CEX=CEX, ...)
    }
    if(type == "summary"){
      oldpar1 <- par(mfrow=c(1,2), pty='s', mar = c(2, 4.25, 0.25, 0.25), 
                     oma = c(0, 0, 0, 0), mgp = c(1, 1, 0), cex.axis=CEX)
      on.exit(par(oldpar1))
      
      plot.psi(out=out, CEX=CEX, ...)
      plot.h.E(out=out, CEX=CEX, ...)
    } 
    if(type == "Qsets"){
      if(missing(col.data)) col.data <- NULL
      if(missing(col.Qfull)) col.Qfull <- NULL
      if(missing(xlim)) xlim <- NULL
      if(missing(ylim)) ylim <- NULL
      
      Plot.Qset.E(obj=out, mar1=mar1, mar2=mar2, dep=dep, est.out=est.out, P=P,   
                  xlim=xlim, ylim=ylim, xlab=labels[1], ylab=labels[2], 
                  col.data=col.data, col.Qfull=col.Qfull, CEX=1.5, ...)
    }
    
  }
  


}


summary_ExtDep <- function(object, mcmc, burn, cred=0.95, plot=FALSE, ...) {
  
  if(missing(object)){
    w <- seq(from=0.0001, to=0.9999, length=100)
  }else{
    w <- object  
  }
  
  nsim <- mcmc$nsim
  ww <- cbind(w, 1-w)
  
  alpha <- (1 - cred)/2
  
  # posterior k
  qk <- quantile(mcmc$k[burn:nsim], c(alpha, 0.5, 1-alpha), type=3)
  k.median <- qk[2]
  k.up <- qk[3]
  k.low <- qk[1]
  
  bpg <- NULL
  for(i in 1:(max(mcmc$k, na.rm=TRUE)+1)) bpg[[i]] <- bp(ww, i)
  
  # h(w) = k sum_j^k-1 (eta_{j+1} - eta_j) bj(w,k-1)
  # A(w) = sum_j^k+1 beta_j bj(w,k+1)
  
  ngrid <- nrow(bpg[[1]]) 
  iters <- nsim-burn+1 
  eta.diff_post <- beta_post <- NULL
  h_post <- A_post <- matrix(NA,iters, ngrid)
  p0_post <- p1_post <- numeric(iters)
  
  if("mar1" %in% names(mcmc)){
    mar1_post <- mar2_post <- matrix(NA, iters, ncol(mcmc$mar1)) # Matrix for the marginal parameters
  }
  
  for(i in burn:nsim){
    l <- i-burn+1
    ki <- mcmc$k[i]
    etai <- mcmc$eta[i,1:(ki+1)]
    eta.diff_post[[l]] <- diff(etai)
    h_post[l,] <- ki * c(bpg[[ki-1]] %*% eta.diff_post[[l]])
    p0_post[l] <- etai[1]
    p1_post[l] <- 1-etai[ki+1]
    beta_post[[l]] <- net(etai, from='H')$beta
    A_post[l,] <- c(bpg[[ki+1]] %*% beta_post[[l]])
    if("mar1" %in% names(mcmc)){ # If the marginal parameters have been fitted
      mar1_post[l,] <- mcmc$mar1[i,]
      mar2_post[l,] <- mcmc$mar2[i,]
    }
  }
  
  # Pointwise credibility bands and mean cuve of h(w)
  h.mean <- apply(h_post, 2, mean)#colMeans(h_post)
  h.up <- apply(h_post, 2, quantile, 1-alpha, type=3)
  h.low <- apply(h_post, 2, quantile, alpha, type=3)
  
  A.mean <- apply(A_post, 2, mean)#colMeans(A_post)
  A.up <- apply(A_post, 2, quantile, 1-alpha, type=3)
  A.low <- apply(A_post, 2, quantile, alpha, type=3)
  
  p0.mean <- mean(p0_post)
  p0.up <- quantile(p0_post, 1-alpha, type=3)
  p0.low <- quantile(p0_post, alpha, type=3)
  
  p1.mean <- mean(p1_post)
  p1.up <- quantile(p1_post, 1-alpha, type=3)
  p1.low <- quantile(p1_post, alpha, type=3)
  
  if("mar1" %in% names(mcmc)){
    mar1.mean <- apply(mar1_post, 2, mean) 
    mar1.up <- apply(mar1_post, 2, quantile, 1-alpha, type=3) 
    mar1.low <- apply(mar1_post, 2, quantile, alpha, type=3) 
    
    mar2.mean <- apply(mar2_post, 2, mean) 
    mar2.up <- apply(mar2_post, 2, quantile, 1-alpha, type=3) 
    mar2.low <- apply(mar2_post, 2, quantile, alpha, type=3)   
  }
  
  if("mar1" %in% names(mcmc)){
    out <- list(k.median=k.median, k.up=k.up, k.low=k.low,
                h.mean=h.mean, h.up=h.up, h.low=h.low,
                A.mean=A.mean, A.up=A.up, A.low=A.low,
                p0.mean=p0.mean, p0.up=p0.up, p0.low=p0.low,
                p1.mean=p1.mean, p1.up=p1.up, p1.low=p1.low,
                mar1.mean=mar1.mean, mar1.up=mar1.up, mar1.low=mar1.low,
                mar2.mean=mar2.mean, mar2.up=mar2.up, mar2.low=mar2.low,
                A_post=A_post, h_post=h_post, 
                eta.diff_post=eta.diff_post, 
                beta_post=beta_post,
                mar1_post=mar1_post, mar2_post=mar2_post,
                p0_post=p0_post, p1_post=p1_post,
                w=w, burn=burn)
  }else{
    out <- list(k.median=k.median, k.up=k.up, k.low=k.low,
                h.mean=h.mean, h.up=h.up, h.low=h.low,
                A.mean=A.mean, A.up=A.up, A.low=A.low,
                p0.mean=p0.mean, p0.up=p0.up, p0.low=p0.low,
                p1.mean=p1.mean, p1.up=p1.up, p1.low=p1.low,
                A_post=A_post, h_post=h_post, 
                eta.diff_post=eta.diff_post, 
                beta_post=beta_post,
                p0_post=p0_post, p1_post=p1_post,
                w=w, burn=burn)
  }
  
  if(plot)
    plot_ExtDep.np(type = "summary", out=mcmc, summary.mcmc=out, burn=burn, ...)
  
  return(out)
  
}

###
### Hidden functions for Bayesian estimation
###

func <- function(y, conf=0.90){
  alpha <- (1-conf)/2
  y <- y[is.finite(y)]
  return(c(quantile(y, alpha), mean(y), quantile(y, 1-alpha)))
}


