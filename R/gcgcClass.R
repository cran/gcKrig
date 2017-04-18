#### In order not to call geoR, just copy this simple function here to be used
"materncopy" <-
  function (u, phi, kappa)
  {
    if(is.vector(u)) names(u) <- NULL
    if(is.matrix(u)) dimnames(u) <- list(NULL, NULL)
    uphi <- u/phi
    uphi <- ifelse(u > 0,
                   (((2^(-(kappa-1)))/ifelse(0, Inf,gamma(kappa))) *
                      (uphi^kappa) *
                      besselK(x=uphi, nu=kappa)), 1)
    uphi[u > 600*phi] <- 0
    return(uphi)
  }


###########################################################################################################
#### Specify a class of marginal functions with the capability of user-defined marginals
#### gc1: parameter values are given by users, for simulation/computing FHUB purposes
###########################################################################################################
#### Continuous Distributions...



gaussian.gc1 <- function(mean = 0, sd = 1){
  q <- function(p) qnorm(p = p, mean = mean, sd = sd)
  ans <- list()
  ans$margvar <- sd^2
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
    integrate(function(x, order)
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rnorm(n = nrep, mean = mean, sd = sd)
  ans$q <- q
  class(ans) <- c("marginal.gc1")
  return(ans)
}




gm.gc1 <- function(shape = 1, rate = 1){
  q <- function(p) qgamma(p = p, shape = shape, rate = rate)
  ans <- list()
  ans$margvar <- shape/(rate^2)
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rgamma(n = nrep, shape = shape, rate = rate)
  ans$q <- q
  class(ans) <- c("marginal.gc1")
  return(ans)
}




beta.gc1 <- function(shape1 = 1, shape2 = 1){
  q <- function(p) qbeta(p = p, shape1 = shape1, shape2 = shape2)
  ans <- list()
  ans$margvar <- shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1))
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rbeta(n = nrep, shape1 = shape1, shape2 = shape2)
  ans$q <- q
  class(ans) <- c("marginal.gc1")
  return(ans)
}





weibull.gc1 <- function(shape = 1, scale = 1){
  q <- function(p) qweibull(p = p, shape = shape, scale = scale)
  ans <- list()
  ans$margvar <- (scale^2)*(base::gamma(1+2/shape) - base::gamma(1+1/shape)^2)
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rweibull(n = nrep, shape = shape, scale = scale)
  ans$q <- q
  class(ans) <- c("marginal.gc1")
  return(ans)
}


poisson.gc1 <- function(lambda = 1){
  q <- function(p) qpois(p = p, lambda = lambda)
  ans <- list()
  ans$margvar <- lambda
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rpois(n = nrep, lambda = lambda)
  ans$q <- q
  class(ans) <- c("marginal.gc1")
  return(ans)
}


negbin.gc1 <- function(mu = 1, od = 1){
  q <- function(p) qnbinom(p = p, size = 1/od, mu = mu)
  ans <- list()
  ans$margvar <- mu*(1+mu*od)
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rnbinom(n = nrep, size = 1/od, mu = mu)
  ans$q <- q
  class(ans) <- c("marginal.gc1")
  return(ans)
}



binomial.gc1 <- function(size = 1, prob = 0.5){
  q <- function(p) qbinom(p = p, size = size, prob = prob)
  ans <- list()
  ans$margvar <- size*prob*(1-prob)
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rbinom(n = nrep, size = size, prob = prob)
  ans$q <- q
  class(ans) <- c("marginal.gc1")
  return(ans)
}





zip.gc1 <- function(mu = 1, od = 1){
  q <- function(p) qpois(pmax( 0, (p-od/(1+od))/(1-od/(1+od)) ), (1+od)*mu)
  ans <- list()
  ans$margvar <- mu*(1+mu*od)
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) ifelse(rbinom(n=nrep, 1, od/(1+od)), 0, rpois(n=nrep, (1+od)*mu))
  ans$q <- q
  class(ans) <- c("marginal.gc1")
  return(ans)
}


####### Correlaiton Functions in the class of gc1 for simulation purpose



matern.gc1 <- function(range = 0, kappa = 0.5, nugget = 0){
  ans <- list()
  ans$S <- function(D) (1-nugget)*materncopy(D, range, kappa) + nugget*diag(NROW(D))
  class(ans) <- c("corr.gc1")
  return(ans)
}




powerexp.gc1 <- function(range = 0, kappa = 1, nugget = 0){
  ans <- list()
  if(kappa > 2)
    stop("'Kappa' must be between 0 and 2")
  ans$S <- function(D) (1-nugget)*exp(-abs( (D/range)^(kappa) )) + nugget*diag(NROW(D))
  class(ans) <- c("corr.gc1")
  return(ans)
}




spherical.gc1 <- function(range = 0, nugget = 0){
  ans <- list()
  ans$S <- function(D) (1-nugget)*(1 - 1.5*D/range + 0.5*(D/range)^3)*(D<=range)+nugget*diag(NROW(D))
  class(ans) <- c("corr.gc1")
  return(ans)
}



###########################################################################################################
#### Specify a class of marginal functions with the potential of user-defined marginals
#### gc2: parameter values are estimated given data, for computing MLE/Prediction
###########################################################################################################


#### Now Try to Create a class of Marginals starting from negative binomial

poisson.gc2 <- function(link = "log"){
  fm <- poisson( substitute( link ) )
  ans <- list()
  ans$start <- function(y, x, effort) {
     mfit <- suppressWarnings(glm.fit(x, y/effort, family = fm))
     mu0 <- fitted(mfit)
     reg0 <- coef(mfit)
     glb <- qnorm(ppois(y-1, lambda = mu0))
     gub <- qnorm(ppois(y, lambda = mu0))
     names(reg0)[1] <- "Intercept"
    return(reg0)
  }
  ans$nod <- 0
  ans$bounds <- function(y, x, pars, effort) {
      M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
      a <- qnorm(ppois( y-1, lambda = M));  b <- qnorm(ppois( y, lambda = M))
    return(list(lower = a, upper = b))
  }
  ans$pdf <- function(y, x, pars, effort){
      M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
      M <- ifelse(M == Inf, .Machine$double.xmax, M)
      pdf <- dpois(y, lambda = M, log = FALSE)
    return(pdf)
  }
  ans$cdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    cdf <- ppois( y, lambda = M)
    return(cdf)
  }
  ans$fm <- fm
  class(ans) <- c("marginal.gc2")
  return(ans)
}


negbin.gc2 <- function(link = "log"){
  fm <- poisson( substitute( link ) )
  ans <- list()
  ans$start <- function(y, x, effort) {
    mfit <- suppressWarnings(glm.fit(x, y/effort, family = fm))
    reg0 <- coef(mfit); mu0 <- fitted(mfit)
    od0 <- max(10*.Machine$double.eps, mean(((y-mu0)^2-mu0)/mu0^2))
    glb <- qnorm(pnbinom(y-1, size = 1/od0, mu = mu0))
    gub <- qnorm(pnbinom(y, size = 1/od0, mu = mu0))
    names( od0 ) <- "overdispersion"
    names(reg0)[1] <- "Intercept"
    return(c(reg0, od0))
  }
  ans$nod <- 1
  ans$bounds <- function(y, x, pars, effort) {
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    S <- 1/pars[ncol(x)+1]
    a <- qnorm(pnbinom( y-1, size = S, mu = M))
    b <- qnorm(pnbinom( y, size = S, mu = M))
    return(list(lower = a, upper = b))
  }
  ans$pdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    S <- 1/pars[ncol(x)+1]
    pdf <- dnbinom(y, size = S, mu = M, log = FALSE)
    return(pdf)
  }
  ans$cdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    S <- 1/pars[ncol(x)+1]
    cdf <- pnbinom(y, size = S, mu = M)
    return(cdf)
  }
  ans$fm <- fm
  class(ans) <- c("marginal.gc2")
 return(ans)
}


binomial.gc2 <- function(link = "logit"){
   fm <- binomial( substitute( link ) )
   ans <- list()
   ans$start <- function(y, x, effort) {
     mfit <- suppressWarnings(glm.fit(x, y/(effort), family = fm))
     reg0 <- coef(mfit)
     glb <- qnorm(pbinom(y-1, size = effort, prob = fitted(mfit)))
     gub <- qnorm(pbinom(y, size = effort, prob = fitted(mfit)))
     names(reg0)[1] <- "Intercept"
     return(reg0)
   }
   ans$nod <- 0
   ans$bounds <- function(y, x, pars, effort) {
     p <- fm$linkinv(pars[1:ncol(x)]%*%t(x))
     a <- qnorm(pbinom( y-1, size = effort, prob = p))
     b <- qnorm(pbinom( y, size = effort, prob = p))
     return(list(lower = a, upper = b))
   }

   ans$pdf <- function(y, x, pars, effort){
     p <- fm$linkinv(pars[1:ncol(x)]%*%t(x))
     pdf <- dbinom(y, size = effort, prob = p, log = FALSE)
     return(pdf)
   }

   ans$cdf <- function(y, x, pars, effort){
     p <- fm$linkinv(pars[1:ncol(x)]%*%t(x))
     cdf <- pbinom(y, size = effort, prob = p)
     return(cdf)
   }
   ans$fm <- fm
   class(ans) <- c("marginal.gc2")
  return(ans)
}


zip.gc2 <- function(link = "log"){
  fm <- poisson( substitute( link ) )
  ans <- list()
  ans$start <- function(y, x, effort) {
    mfit <- suppressWarnings(glm.fit(x, y/effort, family = fm))
    reg0 <- coef(mfit); mu0 <- fitted(mfit)
    od0 <- max(10*.Machine$double.eps, mean(((y-mu0)^2-mu0)/mu0^2))
    glb <- qnorm((y >= 1)*od0/(1+od0)+ppois(y-1, (1+od0)*mu0)/(1+od0))
    gub <- qnorm((y >= 0)*od0/(1+od0)+ppois(y, (1+od0)*mu0)/(1+od0))
    names(reg0)[1] <- "Intercept"
    names(od0) <- "overdispersion"
    return(c(reg0, od0))
  }
  ans$nod <- 1
  ans$bounds <- function(y, x, pars, effort) {
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    od <- pars[ncol(x)+1]
    a <- qnorm( (y>=1)*od/(1+od) + ppois(y-1, lambda = (1+od)*M)/(1+od) )
    b <- qnorm( (y>=0)*od/(1+od) + ppois(y, lambda = (1+od)*M)/(1+od) )
    return(list(lower = a, upper = b))
  }
  ans$pdf <- function(y, x, pars, effort){
   M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
   M <- ifelse(M == Inf, .Machine$double.xmax, M)
   od <- pars[ncol(x)+1]
   pdf <-  (y==0)*od/(1+od) + dpois(y, (1+od)*M)/(1+od)
   return(pdf)
  }
  ans$cdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    od <- pars[ncol(x)+1]
    cdf <- (y >= 0)*od/(1+od)+ppois(y, (1+od)*M)/(1+od)
    return(cdf)
  }
  ans$fm <- fm
  class(ans) <- c("marginal.gc2")
  return(ans)
}



matern.gc2 <- function(kappa = 0.5, nugget = TRUE){
  ans <- list()
    if(nugget == TRUE){
      ans$nug <- 1
      ans$npar.cor <- 2
      ans$start <- function(D) {
        corstart <- c(median(D)/2,0.2)
        names(corstart) <- c("range","nugget")
        return(corstart)
      }
      ans$corr <- function(corrpar,D){
        S <- (1 - corrpar[2]) * materncopy(D, corrpar[1], kappa) + corrpar[2]*diag(NROW(D))
        return(S)
      }
    }else{
      ans$nug <- 0
      ans$npar.cor <- 1
      ans$start <- function(D)  {
        corstart <- median(D)/2
        names(corstart) <- "range"
        corstart
      }
      ans$corr <- function(corrpar,D){
        S <- materncopy(D, corrpar, kappa)
        return(S)
      }
    }
  class(ans) <- c("corr.gc2")
  return(ans)
}



powerexp.gc2 <- function(kappa = 1, nugget = TRUE){
  ans <- list()
  if(kappa > 2)
    stop("'Kappa' must be between 0 and 2")
  if(nugget == TRUE){
    ans$nug <- 1
    ans$npar.cor <- 2
    ans$start <- function(D) {
      corstart <- c(median(D)/2, 0.2)
      names(corstart) <- c("range","nugget")
      return(corstart)
    }
    ans$corr <- function(corrpar, D ){
      S <- (1 - corrpar[2])*exp(-abs( (D/corrpar[1])^(kappa))) + corrpar[2]*diag(NROW(D))
      return(S)
    }
  }else{
    ans$nug <- 0
    ans$npar.cor <- 1
    ans$start <- function(D) {
      corstart <- median(D)/2
      names(corstart) <- "range"
      return(corstart)
    }
    ans$corr <- function(corrpar, D ){
      S <- exp(-abs((D/corrpar[1])^(kappa)))
      return(S)
    }
  }
  class(ans) <- c("corr.gc2")
  return(ans)
}


spherical.gc2 <- function(nugget = TRUE){
  ans <- list()
  if(nugget == TRUE){
    ans$nug <- 1
    ans$npar.cor <- 2
    ans$start <- function(D) {
      corstart <- c(median(D)/2, 0.2)
      names(corstart) <- c("range","nugget")
      return(corstart)
    }
    ans$corr <- function(corrpar, D ){
      S <- (1 - corrpar[2])*(1 - 1.5*D/corrpar[1] +
                               0.5*(D/corrpar[1])^3)*(D<=corrpar[1]) + corrpar[2]*diag(NROW(D))
      return(S)
    }
  }else{
    ans$nug <- 0
    ans$npar.cor <- 1
    ans$start <- function(D) {
      corstart <- median(D)/2
      names(corstart) <- "range"
      return(corstart)
    }
    ans$corr <- function(corrpar, D ){
      S <- (1 - 1.5*D/corrpar[1] + 0.5*(D/corrpar[1])^3)*(D<=corrpar[1])
      return(S)
    }
  }
  class(ans) <- c("corr.gc2")
  return(ans)
}


