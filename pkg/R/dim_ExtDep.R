dim_ExtDep <- function(model, par=NULL, dim=NULL){
  
  models <- c("PB", "HR", "ET", "EST", "TD", "AL")
  if(!any(model ==  models)){ stop("model wrongly specified")}
  
  if(is.null(par) && is.null(dim)){ stop("par or dim needs to be specified")}
  
  if(!is.null(par) && is.null(dim)){
    if(!is.vector(par)){ stop("par is not a vector")}
    
    if(model == "HR"){
      dim <- round(uniroot( function(y) choose(y,2) -length(par), interval=c(1,50)  )$root)
      if(length(par) != (choose(dim,2)) ){stop("Incorrect parameter dimension for model specified")}
    }else if(model == "PB" || model == "ET"){
      dim <- round(uniroot( function(y) choose(y,2) + 1 -length(par), interval=c(1,50)  )$root)
      if(length(par) != (choose(dim,2)+1) ){stop("Incorrect parameter dimension for model specified")}
    }else if(model == "EST"){
      dim <- round(uniroot( function(y) choose(y,2) + y + 1 -length(par), interval=c(1,50)  )$root)
      if(length(par) != (choose(dim,2)+dim+1) ){stop("Incorrect parameter dimension for model specified")}
    }else if(model == "TD"){
      dim <- length(par)
    }else if(model ==  "AL"){
      dim <- round(uniroot( function(y) 2^(y-1) * (y+2) - (2*y+1) -length(par), interval=c(1,50)  )$root)
      if(length(par) != (2^(dim-1)*(dim+2) - (2*dim+1) ) ){stop("Incorrect parameter dimension for model specified")}
    }
    
    return(dim)
  }
  
  if(is.null(par) && !is.null(dim)){
    if(model == "HR"){
      length.par = choose(dim,2)
    }else if(model == "PB" || model == "ET"){
      length.par =  choose(dim,2) + 1
    }else if(model == "EST"){
      length.par =  choose(dim,2) + dim + 1
    }else if(model == "TD"){
      length.par = dim
    }else if(model ==  "AL"){
      length.par = 2^(dim-1)*(dim+2) - (2*dim+1) 
    }
    
    return(length.par)
  }
  
  if(!is.null(par) && !is.null(dim)){
    if(model == "HR"){
      length(par) == choose(dim,2)
    }else if(model == "PB" || model == "ET"){
      length(par) ==  choose(dim,2) + 1
    }else if(model == "EST"){
      length(par) ==  choose(dim,2) + dim + 1
    }else if(model == "TD"){
      length(par) == dim
    }else if(model ==  "AL"){
      length(par) == 2^(dim-1)*(dim+2) - (2*dim+1) 
    }
    
  }
  
}
