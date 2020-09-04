as.mcmc.blca.collapsed <-
  function(x, ...){
    
    if(inherits(x, "blca.collapsed", TRUE) == 0) stop("Invalid Class, must be blca.collapsed object")
    else{
      y <- x$samples$G
      as.mcmc(y, ...)
    }
  }
