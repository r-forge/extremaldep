#######################################################
### Authors: Giulia Marcon and Simone Padoan        ###
### Emails: giulia.marcon@phd.unibocconi.it,        ###
### simone.padoan@unibocconi.it                     ###
### Institution: Department of Decision Sciences,   ###
### University Bocconi of Milan                     ###
### File name: plots.r                              ###
### Description:                                    ###
### This file provides the function for graphical   ###
### representations of the Bayesian Inference for   ###
### the Extremal dependence, in Bernstein form      ###
### (see Marcon et al., 2016)                       ###
### Last change: 15/08/2016                         ###
#######################################################


plot.bbeed <- function(x, type = c("summary", "returns", "A", "h", "pm", "k"), 
                       mcmc, summary.mcmc, nsim, burn, y, probs, CEX=1.5, A_true, h_true,
                       labels=c(expression(y[1]),expression(y[2])), ...){
  
###################### START PLOT RETURNS ####################

################################################################################
# INPUT:                                                                     ###
# y is a (m x 2)-dimensional matrix of the thresholds                        ###
# probs is the output obtained with returns function, i.e. vector of       ###
#          returns values                                                    ###
################################################################################ 

  plot.returns <- function(y, probs, CEX=1.5,
                         labels=c(expression(y[1]),expression(y[2])), ...){
  
  op1 <- par(mar = c(1, 1, 0, 0), oma = c(3, 4, 0.5, 0.5))
  on.exit(par(op1))
  op2 <- par(mgp = c(1, 1, 0),cex.axis=CEX)
  on.exit(par(op2))
  
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
  
  mtext('w',side=1,line=3,cex=CEX)
  mtext('h(w)',side=2,line=3,cex=CEX)  
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
  mtext('k',side=1,line=3,cex=CEX)
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

if(type == "returns"){
  plot.returns(y=y, probs=probs, CEX=CEX, labels=labels, ...)
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
  PriorVSPosterior.pm(mcmc=mcmc, nsim=nsim, burn=burn, CEX=CEX, ...)
}
if(type == "k"){
  PriorVSPosterior.k(mcmc=mcmc, nsim=nsim, burn=burn, CEX=CEX, ...)
}
if(type == "summary"){
  oldpar1 <- par(mfrow=c(2,2), pty='s')
  on.exit(par(oldpar1))
  oldpar2 <- par(mar = c(4, 0, 0.5, 0), oma = c(3, 4, 0.5, 0.5))
  on.exit(par(oldpar2))
  oldpar3 <- par(mgp = c(1, 1, 0), cex.axis=CEX)
  on.exit(par(oldpar3))
  
  plot.A(w=x, summary.mcmc=summary.mcmc, ...)
  PriorVSPosterior.k(mcmc=mcmc, nsim=nsim, burn=burn, ...)
  plot.h(w=x, summary.mcmc=summary.mcmc, ...)
  PriorVSPosterior.pm(mcmc=mcmc, nsim=nsim, burn=burn, ...)
}

}
###################### END PLOT PRIOR VS POSTERIOR p0 ####################
