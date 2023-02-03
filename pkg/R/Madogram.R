#######################################################
### Authors: Giulia Marcon and Simone Padoan        ###
### Emails: giulia.marcon@phd.unibocconi.it,        ###
### simone.padoan@unibocconi.it                     ###
### Institution: Department of Decision Sciences,   ###
### University Bocconi of Milan                     ###
### File name: Madogram.r                           ###
### Description:                                    ###
### This file enables to compute the multivariate   ###
### madogram as proposed in Marcon et al. (2014)    ###
### Last change: 14/07/2015                         ###
#######################################################

# The routine estimates the Pickands function using
# the multivariate madogram
# Multivariate case

madogram <- function (w, data, margin = c("emp", "est", "exp", "frechet", "gumbel")) 
{
  lmadogram <- function(w, data, margin) {
    sumdata <- dim(data)
    d <- sumdata[2]
    if (!is.matrix(w)) 
      w <- matrix(w, ncol = d)
    sumw <- dim(w)
    if (sumw[2] != d) 
      stop("`x' must be a vector/matrix with `d' elements/columns")
    if( length(margin)>1 ) 
      margin = "emp"
    ans <- numeric(sumw[1])
    if (margin == "emp") {
      data_emp <- apply(data, 2, rank, na.last = "keep")
      nasm <- apply(data_emp, 2, function(x) sum(!is.na(x)))
      data_emp <- data_emp/rep(nasm , each = nrow(data_emp))
      Fdata <- data_emp
    }
    if(margin=="est"){
      par <- NULL
      Fdata <- data
      for(i in 1:d){
        par[[i]] <- fGEV(data[,i])
        Fdata[,i] <- pgev(data[,i], loc=par[[i]]$param[1], 
                          scale=par[[i]]$param[2], shape=par[[i]]$param[3])
      }
    }
    if (margin == "exp") {
      Fdata <- apply(data, 2, pexp)
    }
    if (margin == "frechet") {
      Fdata <- apply(data, 2, pfrechet)
    }
    if (margin == "gumbel") {
      Fdata <- apply(data, 2, pgumbel)
    }
    powerw <- function(j, xx, w, d) 
      sapply(c(1:d), function(i, x, w) x[, i]^(1/w[, i]), xx, t(w[j, ]))
    u <- lapply(c(1:sumw[1]), powerw, Fdata, w, d)
    ma <- sapply(c(1:sumw[1]), function(i, u) apply(u[[i]], 1, max), u)
    me <- sapply(c(1:sumw[1]), function(i, u) rowMeans(u[[i]]), u)
    mm <- ma - me
    ans <- colMeans(mm)
    return(ans)
  }
  ans <- lmadogram(w, data, margin)
  sumdata <- dim(data)
  d <- sumdata[2]
  if (!is.matrix(w)) 
    w <- matrix(w, ncol = d)
  W <- w/(1 + w)
  cst <- rowMeans(W)
  A <- ((cst + ans)/(1 - cst - ans))
  return(A)
}