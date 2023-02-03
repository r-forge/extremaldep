rExtDepSpat <- function(n, coord, model="SCH", cov.mod = "whitmat", grid = FALSE,
                        control = list(), cholsky = TRUE, ...){

  
    models <- c("SMI","GG", "BR", "SCH", "ET", "EST")
    if(!any(model ==  models)){ stop("model wrongly specified")}
    
    if(model == "SMI"){
      cov.mod <- "gauss"  
    }else if(model == "BR"){
      cov.mond <- "brown"
    }else if(model %in% c("GG", "SCH", "ET", "EST")){
      if(!(cov.mod %in% c("whitmat","cauchy","powexp","bessel"))){
        stop("'cov.mod' must be one of 'whitmat', 'cauchy', 'powexp' or 'bessel'.")
      }
    }
    
    # if (!(cov.mod %in% c("whitmat","cauchy","powexp","bessel",
    #                      "iwhitmat", "icauchy", "ipowexp", "ibessel",
    #                      "gwhitmat", "gcauchy", "gpowexp", "gbessel",
    #                      "twhitmat", "tcauchy", "tpowexp", "tbessel",
    #                      "swhitmat", "scauchy", "spowexp", "sbessel",
    #                      "brown")))
    #     stop("'cov.mod' must be one of '(i/g/t/s)whitmat', '(i/g/t/s)cauchy', '(i/g/t/s)powexp', '(i/g/t/s)bessel' or 'brown'")

    if (!is.null(control$method) && !(control$method %in% c("direct", "tbm", "circ", "exact")))
        stop("the argument 'method' for 'control' must be one of 'exact', 'direct', 'tbm' and 'circ'")

    # if (cov.mod == "brown")
    #     model <- "Brown-Resnick"
    # 
    # else if (cov.mod %in% c("whitmat","cauchy","powexp","bessel"))
    #     model <- "Schlather"
    # 
    # else {
    #     if (substr(cov.mod, 1, 1) == "i")
    #         model <- "iSchlather"
    # 
    #     else if (substr(cov.mod, 1, 1) == "g")
    #         model <- "Geometric"
    #         
    #     else if (substr(cov.mod, 1, 1) == "s")
    #         model <- "Extremal-Skew-t"
    # 
    #     else
    #         model <- "Extremal-t"
    # 
    #     cov.mod <- substr(cov.mod, 2, 10)
    # }

    dist.dim <- ncol(coord)

    if (is.null(dist.dim))
        dist.dim <- 1

    if (dist.dim > 2)
        stop("Currently this function is only available for R or R^2")

    if ((dist.dim == 1) && grid){
        warning("You cannot use 'grid = TRUE' in dimension 1. Ignored.")
        grid <- FALSE
    }

    if (model == "SCH"){
        if (!all(c("nugget", "range", "smooth") %in% names(list(...))))
            stop("You must specify 'nugget', 'range', 'smooth'")

        nugget <- list(...)$nugget
        range <- list(...)$range
        smooth <- list(...)$smooth
    }

    else if (model == "GG") {
        if (!all(c("sigma2", "nugget", "range", "smooth") %in% names(list(...))))
            stop("You must specify 'sigma2', 'nugget', 'range', 'smooth'")

        nugget <- list(...)$nugget
        range <- list(...)$range
        smooth <- list(...)$smooth
        sigma2 <- list(...)$sigma2
    }

    else if (model == "ET"){
        if (!all(c("DoF", "nugget", "range", "smooth") %in% names(list(...))))
            stop("You must specify 'DoF', 'nugget', 'range', 'smooth'")

        nugget <- list(...)$nugget
        range <- list(...)$range
        smooth <- list(...)$smooth
        DoF <- list(...)$DoF
    }
    
    else if (model == "EST"){
        if (!all(c("DoF", "nugget", "range", "smooth", "alpha", "acov1", "acov2") %in% names(list(...))))
            stop("You must specify 'DoF', 'nugget', 'range', 'smooth', 'alpha','acov1', 'acov2' ")

        nugget <- list(...)$nugget
        range <- list(...)$range
        smooth <- list(...)$smooth
        DoF <- list(...)$DoF
        alpha <- list(...)$alpha
        acov1 <- list(...)$acov1
        acov2 <- list(...)$acov2        
    }

    else if (model == "BR"){
        range <- list(...)$range
        smooth <- list(...)$smooth
    }

    if (dist.dim !=1){
        n.site <- nrow(coord)
        coord.range <- apply(coord, 2, range)
        center <- colMeans(coord.range)
        edge <- max(apply(coord.range, 2, diff))
    }

    else {
        n.site <- length(coord)
        coord.range <- range(coord)
        center <- mean(coord.range)
        edge <- diff(coord.range)
    }


    cov.mod <- switch(cov.mod, "gauss" = "gauss", "whitmat" = 1, "cauchy" = 2,
                      "powexp" = 3, "bessel" = 4)

    if (grid){
        n.eff.site <- n.site^dist.dim
        reg.grid <- .isregulargrid(coord[,1], coord[,2])
        steps <- reg.grid$steps
        reg.grid <- reg.grid$reg.grid
    }

    else
        n.eff.site <- n.site

    ans <- rep(-1e10, n * n.eff.site)
    ans2 <- rep(0L, n * n.eff.site)

    ##Identify which simulation technique is the most adapted or use the
    ##one specified by the user --- this is useless for the Smith model.
    if (is.null(control$method)){
        if (n.eff.site < 500)## <<-- If fixed it to 500 but need to modify it later maybe !!!
            method <- "exact"

        else if (grid && reg.grid)
            method <- "circ"

        else if ((length(ans) / n) > 600)
            method <- "tbm"

        else
            method <- "direct"
    } else {
        method <- control$method
    }


    if (method == "tbm"){
        if (is.null(control$nlines))
            nlines <- 1000

        else
            nlines <- control$nlines
    }

    if (model == "SCH"){

        if (is.null(control$uBound))
            uBound <- 3.5

        else
            uBound <- control$uBound

        if (method == "direct")
            ans <- .C("rschlatherdirect", as.double(coord), as.integer(n), as.integer(n.site),
                      as.integer(dist.dim), as.integer(cov.mod), grid, as.double(nugget),
                      as.double(range), as.double(smooth), as.double(uBound), ans = ans, ans2 = ans2, NAOK=TRUE)[c("ans","ans2")]

        else if (method == "exact")
            ans <- .C("rschlatherexact", as.double(coord), as.integer(n), as.integer(n.site), as.integer(dist.dim),
                      as.integer(cov.mod), as.integer(grid), as.double(nugget), as.double(range), as.double(smooth),
                      ans = ans, ans2 = ans2, NAOK=TRUE)[c("ans","ans2")]

        else if (method == "circ")
            ans <- .C("rschlathercirc", as.integer(n), as.integer(n.site), as.double(steps),
                      as.integer(dist.dim), as.integer(cov.mod), as.double(nugget), as.double(range),
                      as.double(smooth),  as.double(uBound), ans = ans, NAOK=TRUE)$ans

        else
            ans <- .C("rschlathertbm", as.double(coord), as.integer(n), as.integer(n.site),
                      as.integer(dist.dim), as.integer(cov.mod), grid, as.double(nugget),
                      as.double(range), as.double(smooth), as.double(uBound), as.integer(nlines),
                      ans = ans, NAOK=TRUE)$ans
        
    }

    else if (model == "GG"){

        if (is.null(control$uBound))
            uBound <- exp(3.5 * sqrt(sigma2) - 0.5 * sigma2)

        else
            uBound <- control$uBound

        if (method == "direct")
            ans <- .C("rgeomdirect", as.double(coord), as.integer(n), as.integer(n.site),
                      as.integer(dist.dim), as.integer(cov.mod), grid, as.double(sigma2),
                      as.double(nugget), as.double(range), as.double(smooth),
                      as.double(uBound), ans = ans, ans2 = ans2, NAOK=TRUE)[c("ans","ans2")]
                      
        else if (method == "exact")
            ans <- .C("rgeomexact", as.double(coord), as.integer(n), as.integer(n.site), as.integer(dist.dim),
                      as.integer(cov.mod), as.integer(grid), as.double(sigma2), as.double(nugget), as.double(range), 
                      as.double(smooth), ans = ans, ans2 = ans2, NAOK=TRUE)[c("ans","ans2")]

        else if (method == "circ")
            ans <- .C("rgeomcirc", as.integer(n), as.integer(n.site), as.double(steps),
                      as.integer(dist.dim), as.integer(cov.mod), as.double(sigma2),
                      as.double(nugget), as.double(range), as.double(smooth),  as.double(uBound),
                      ans = ans, NAOK=TRUE)$ans

        else
            ans <- .C("rgeomtbm", as.double(coord), as.integer(n), as.integer(n.site),
                      as.integer(dist.dim), as.integer(cov.mod), grid, as.double(sigma2),
                      as.double(nugget), as.double(range), as.double(smooth), as.double(uBound),
                      as.integer(nlines), ans = ans, NAOK=TRUE)$ans
    }

    else if (model == "ET"){
        if (is.null(control$uBound))
            uBound <- 3^DoF

        else
            uBound <- control$uBound

        if (method == "direct")
            ans <- .C("rextremaltdirect", as.double(coord), as.integer(n), as.integer(n.site),
                      as.integer(dist.dim), as.integer(cov.mod), grid, as.double(nugget),
                      as.double(range), as.double(smooth), as.double(DoF), as.double(uBound),
                      ans = ans, ans2 = ans2, NAOK=TRUE)[c("ans","ans2")]

        else if (method == "exact")
            ans <- .C("rextremaltexact", as.double(coord), as.integer(n), as.integer(n.site),
                      as.integer(dist.dim), as.integer(cov.mod), grid, as.double(nugget),
                      as.double(range), as.double(smooth), as.double(DoF), as.integer(cholsky), 
                      ans = ans, ans2 = ans2, NAOK=TRUE)[c("ans","ans2")]


        else if (method == "circ")
            ans <- .C("rextremaltcirc", as.integer(n), as.integer(n.site), as.double(steps),
                      as.integer(dist.dim), as.integer(cov.mod), as.double(nugget), as.double(range),
                      as.double(smooth), as.double(DoF), as.double(uBound), ans = ans, NAOK=TRUE)$ans

        else
            ans <- .C("rextremalttbm", as.double(coord), as.integer(n), as.integer(n.site),
                      as.integer(dist.dim), as.integer(cov.mod), grid, as.double(nugget),
                      as.double(range), as.double(smooth), as.double(DoF), as.double(uBound),
                      as.integer(nlines), ans = ans, NAOK=TRUE)$ans
    }
    
    else if (model == "EST"){
        if (is.null(control$uBound))
            uBound <- 3^DoF

        else
            uBound <- control$uBound

        if(length(alpha) != 3) stop("alpha must be of length 3")
        if(length(acov1) != n.eff.site) stop("acov1 has incorrect length")
        if(length(acov2) != n.eff.site) stop("acov2 has incorrect length")
        
        Alpha <- alpha[1] + alpha[2] * acov1 + alpha[3] * acov2
        
        if (method == "direct")
            ans <- .C("rextremalskewtdirect", as.double(coord), as.integer(n), as.integer(n.site),
                      as.integer(dist.dim), as.integer(cov.mod), grid, as.double(nugget),
                      as.double(range), as.double(smooth), as.double(DoF), as.double(Alpha), as.double(uBound),
                      ans = ans, ans2 = ans2, NAOK=TRUE)[c("ans","ans2")]

        else if (method == "exact")
            ans <- .C("rextremalskewtexact", as.double(coord), as.integer(n), as.integer(n.site),
                      as.integer(dist.dim), as.integer(cov.mod), grid, as.double(nugget),
                      as.double(range), as.double(smooth), as.double(DoF), as.double(Alpha), 
                      as.integer(cholsky), ans = ans, ans2 = ans2, NAOK=TRUE)[c("ans","ans2")]

        else stop("Extremal Skew-t only has direct and exact methods")
    }

    else if (model == "BR"){
        coord <- scale(coord, scale = FALSE) + 10^-6## to avoid having the origin

        if (is.null(control$max.sim))
            max.sim <- 1000

        else
            max.sim <- control$max.sim

        if (is.null(control$uBound))
            uBound <- 10

        else
            uBound <- control$uBound

        if (is.null(control$sim.type))
            sim.type <- 1

        else
            sim.type <- control$sim.type

        if (dist.dim == 1)
            bounds <- range(coord)

        else
            bounds <- apply(coord, 2, range)

        if (is.null(control$nPP))
            nPP <- 15

        else
            nPP <- control$nPP

        if (sim.type == 6){
            idx.sub.orig <- getsubregions(coord, bounds, range, smooth, dist.dim)
            n.sub.orig <- length(idx.sub.orig)
        }

        else
            idx.sub.orig <- n.sub.orig <- 0

        if (method == "direct")
            ans <- .C("rbrowndirect", as.double(coord), as.double(bounds),
                      as.integer(n), as.integer(n.site), as.integer(dist.dim),
                      as.integer(grid), as.double(range), as.double(smooth),
                      as.double(uBound), as.integer(sim.type), as.integer(max.sim),
                      as.integer(nPP), as.integer(idx.sub.orig), as.integer(n.sub.orig),
                      ans = ans, ans2 = ans2, NAOK=TRUE)[c("ans","ans2")]

        else if (method == "exact")
            ans <- .C("rbrownexact", as.double(coord), as.integer(n), as.integer(n.site),
                      as.integer(dist.dim), as.integer(grid), as.double(range), as.double(smooth),
                      ans = ans, ans2 = ans2, NAOK=TRUE)[c("ans","ans2")]
    }

    if(is.list(ans)) {
        ans2 <- ans$ans2; ans <- ans$ans
    } else ans2 <- NA
        
    if (grid){
        if (n == 1) {
            ans <- matrix(ans, n.site, n.site)
            ans2 <- matrix(ans2, n.site, n.site)
        }
        else {
            ans <- array(ans, c(n.site, n.site, n))
            ans2 <- array(ans2, c(n.site, n.site, n))
        }
    }

    else {
        ans <- matrix(ans, nrow = n, ncol = n.site)
        ans2 <- matrix(ans2, nrow = n, ncol = n.site)
    }

    return(list(vals = ans, hits = ans2))
}

getsubregions <- function(coord, bounds, range, smooth, dist.dim){

    if (dist.dim == 1){
        coord <- matrix(coord, ncol = 1)
        bounds <- matrix(bounds, ncol = 1)
    }

    h.star <- 2^(1/smooth) * range
    n.windows <- 1 + floor(diff(bounds) / h.star)

    sub.bounds <- list()
    for (i in 1:dist.dim)
        sub.bounds <- c(sub.bounds,
                        list(seq(bounds[1,i], bounds[2, i], length = n.windows[i] + 1)))

    sub.centers <- list()
    for (i in 1:dist.dim)
        sub.centers <- c(sub.centers, list(0.5 * (sub.bounds[[i]][-1] +
                                                  sub.bounds[[i]][-(n.windows + 1)])))

    sub.origins <- as.matrix(expand.grid(sub.centers))

    n.sub.origins <- prod(n.windows)
    idx.sub.orig <- rep(NA, n.sub.origins)
    for (i in 1:n.sub.origins){
        dummy <- sqrt(colSums((t(coord) - sub.origins[i,])^2))
        idx.sub.orig[i] <- which.min(dummy)
    }


    return(idx.sub.orig = idx.sub.orig)
}

###############################################################################
###############################################################################
## Hidden functions
###############################################################################
###############################################################################

.isregulargrid <- function(x, y, tol.x = 1e-4, tol.y = 1e-4){
  ##This function check if the grid defined by x and y is regular
  ##i.e. the spacings along the axis is constant -- but not
  ##necessarily the same for the x and y-axis
  
  x.diff <- diff(x)
  y.diff <- diff(y)
  
  eps.x <- diff(range(x.diff))
  eps.y <- diff(range(y.diff))
  
  if ((eps.x <= tol.x) && (eps.y <= tol.y)){
    reg.grid <- TRUE
    steps <- c(mean(x.diff), mean(y.diff))
  }
  
  else{
    reg.grid <- FALSE
    steps <- NA
  }
  
  return(list(reg.grid = reg.grid, steps = steps))
}

.isregulargrid <- function(x, y, tol.x = 1e-4, tol.y = 1e-4){
  ##This function check if the grid defined by x and y is regular
  ##i.e. the spacings along the axis is constant -- but not
  ##necessarily the same for the x and y-axis
  
  x.diff <- diff(x)
  y.diff <- diff(y)
  
  eps.x <- diff(range(x.diff))
  eps.y <- diff(range(y.diff))
  
  if ((eps.x <= tol.x) && (eps.y <= tol.y)){
    reg.grid <- TRUE
    steps <- c(mean(x.diff), mean(y.diff))
  }
  
  else{
    reg.grid <- FALSE
    steps <- NA
  }
  
  return(list(reg.grid = reg.grid, steps = steps))
}

