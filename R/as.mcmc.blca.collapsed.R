as.mcmc.blca.collapsed <-
  function(x, ...){
    
    if(inherits(x, "blca.collapsed", TRUE) == 0) stop("Invalid Class, must be blca.collapsed object")
    else{
      
      y <- cbind(x$samples$logpost[-(1:x$burn.in)], x$samples$Giter, rowSums(x$samples$var.ind))
      colnames(y)<- c("log-posterior", "number-of-groups", "number-of-variables")
      as.mcmc(y, ...)
    }
  }
