blca.boot <- function( X, G, ncat=NULL, alpha=1,beta=1, delta=1, start.vals= c("single","across"), counts.n=NULL, model.indicator=NULL, fit=NULL, iter=2000,  B=100, verbose=TRUE, verbose.update=100, conv=1e-06, small=1e-10, MAP=TRUE)
  {
    # check if data is simulated 
    if( class(X) == "blca.rand" & !is.matrix(X) ) X <- X$X
    
    #list of returns
    args.passed <- as.list( environment() )
    x <- list()
    x$call <- match.call()
    x$args <- args.passed
    
    # blca.check.missing ...
    miss <- blca.check.missing( X )
    if( miss$missing ){
      warning("Missing values encountered in X: rows with NA have been removed", call.=FALSE)
      X <- na.omit(X)
    } 
    
    # convert X into a numeric matrix and check inputs
    D <- blca.check.data( X, counts.n, ncat )
    
    X <- D$X
    ncat <- D$ncat
    x$args$ncat <- ncat
    x$G <- G
    
    N<-nrow(X) 
    M<-ncol(X) 
    
    if( is.null(model.indicator) )
    {
      model.indicator <- rep(1,M)
    }else if( length(model.indicator) != M ){
      stop("model.indicator must have length ncol(X)")
    }
    
    M.in <- sum( model.indicator ) 
    if( prod( ncat[model.indicator==1] ) <= (M.in+1)*G) stop("maximum numer of classes that should be run for this data is ", floor(prod(ncat[model.indicator==1])/(M.in+1)) )
    
    if(is.null(fit)){
      if(verbose==TRUE) cat("Object 'fit' not supplied: obtaining starting values via blca.em...\n")
      if(is.null(G)) warning("Number of groups must be specified", call.=FALSE )
      xx<-blca.em(X, G, ncat=ncat, iter=1000, alpha=alpha, beta=beta, delta=delta , start.vals= start.vals, model.indicator=model.indicator, for.boot=TRUE, verbose=FALSE, conv=conv) 
      #conv<-xx$eps
      if(verbose==TRUE) cat("Starting values obtained...\n")
    }else{
      if( is.null( fit$for.boot ) ) stop("rerun blca.em with for.boot set to TRUE")
      xx<- fit
      G<- length(xx$classprob)
      #conv<- xx$eps
      if( xx$MAP )
      {
        alpha<- xx$prior$alpha
        beta<- xx$prior$alpha
        delta<- xx$prior$delta
      }
    }
    
    Z.ref <- xx$Z 
    
    out.prior <- blca.check.prior( alpha, beta, delta, G, M, ncat )
    prior.init.type <- out.prior$prior.init.type
    gamma <- out.prior$gamma
    delta <- out.prior$delta	
    
    hparam <- c( delta[1], beta[1] )
    group.probs <- numeric(B*G*N)
    weights <- numeric(B*G)
    variable.probs <- numeric(B*G*sum(ncat))
    log.post <- log.like <- numeric( B )
    
    # indexes for bootstrap samples
    boot.samp.idx <- matrix( sample( 0:(N-1), size=B*N, replace=T ), nrow=N )
    boot.samp.idx <- as.vector( apply( boot.samp.idx, 2, sort ) )
    
    if( is.null(model.indicator) ) model.indicator <- rep(1,M)
    
    if( verbose ) cat("Beginning bootstrapping run...")
    x.vec <- as.vector(t(X))
    
    w <- .C( 		"BLCA_BOOT_FIT", 						as.integer(x.vec),
               as.integer(N),	
               as.integer(M),							as.integer(ncat),
               as.double(hparam),					as.integer(prior.init.type),
               as.double(delta),						as.double(gamma),
               as.integer(B),							as.integer(G),
               as.integer(boot.samp.idx),
               as.double(xx$boot.init$weights),	as.double(xx$boot.init$var.probs),
               group.probs = as.double(group.probs),
               weights = as.double(weights),		variable.probs = as.double(variable.probs),
               logpost = as.double(log.post),
               loglike = as.double(log.like),
               as.integer(model.indicator),		as.integer(MAP),
               as.double(xx$conv),					as.integer(iter),
               as.integer(verbose), 				as.integer(verbose.update), 
               PACKAGE = "BayesLCA"	)
    
    
    if( verbose ) cat("\nFinished Bootstrap run...")
    
    #split these into a list as it will be easier to rearrange from label processing
    var.probs.l <- vector( M.in, mode="list" )
    
    l <- 1
    for( j in 1:M )
    {
      if(j == 1){
        gap <- 0
      }else{
        gap <- B * G * sum( ncat[1:(j-1)] )
      }
      #variable probabilities stacked by group and then iteration
      if( model.indicator[j] )
      {
        var.probs.l[[l]] <- matrix( w$variable.probs[(gap+1):(gap + B*G*ncat[j])] , nrow = B * G, ncol=ncat[j], byrow=TRUE )
        l <- l+1
      }
    }	
    
    weights.mat <- matrix( w$weights , nrow=B, ncol=G, byrow=TRUE )	
    
    # correct for potential label swapping
    Zarr <- array( w$group.probs, dim=c(N,G,B) )
    
    if( verbose ) cat("\nPost processing bootstrap output...\n")
    
    boot.samp.idx <- matrix( boot.samp.idx, nrow=N ) + 1	
    for( b in 1:B )
    {
      rl <- matchClasses( t(Z.ref[ boot.samp.idx[,b], ]) %*% Zarr[,,b], method="exact", verbose=FALSE )
      weights.mat[ b, ] <- weights.mat[ b, rl ]
      l <- 1
      for( j in 1:M  )
      {
        it <- (b-1)*G
        if( model.indicator[j] )
        {
          var.probs.l[[l]][ it + 1:G , ] = var.probs.l[[l]][ it + rl , ]
          l <- l+1
        }
      }
    }
    
    # get weights to order in size of cluster
    x$classprob <- apply( weights.mat, 2, mean )	
    o <- order(x$classprob, decreasing=TRUE)
    x$classprob <- x$classprob[o]
    
    x$classprob.sd <- apply( weights.mat, 2, sd )
    x$classprob.sd <- x$classprob.sd[o]
    
    # compute for results
    v.probs <- list( mean=vector(M.in, mode="list"), sd=vector(M.in, mode="list") )
    tt <- seq( 0, G*(B-1), by=G ) # corrected this to be G*(B-1)
    l <- 1
    for( j in 1:M )
    {
      if( model.indicator[j] )
      {
        v.probs$mean[[l]]  <- matrix( nrow = G, ncol = ncat[j] ) 
        v.probs$sd[[l]] <- matrix( nrow = G, ncol = ncat[j] )
        for( g in 1:G )
        {
          bm <- var.probs.l[[l]][ tt+g , ]
          v.probs$mean[[l]][g,] <- apply( bm, 2, mean )
          v.probs$sd[[l]][g,] <- apply( bm, 2, sd )
        }
        v.probs$mean[[l]] <- v.probs$mean[[l]][o,]
        v.probs$sd[[l]] <- v.probs$sd[[l]][o,]
        rownames( v.probs$mean[[l]] ) <- paste( "Group", 1:G )
        colnames( v.probs$mean[[l]] ) <- paste("Cat",0:(ncat[j]-1) )
        rownames( v.probs$sd[[l]] ) <- paste( "Group", 1:G )
        colnames( v.probs$sd[[l]] ) <- paste("Cat",0:(ncat[j]-1) )
        l <- l+1
      }
    }		
    
    if( is.null(colnames(X)) )
    {
      names( v.probs$mean ) <- names( v.probs$sd ) <- paste("Variable", which( model.indicator==1 ) )
    }else{
      names( v.probs$mean ) <- names( v.probs$sd ) <- colnames(X)[  which( model.indicator == 1 ) ]
    }	
    
    #compile the list to return
    x$itemprob <- v.probs$mean
    x$itemprob.sd <- v.probs$sd
    
    vec.itemprobs <- unlist(x$itemprob)
    vec.itemprobs.sd <- unlist(x$itemprob.sd)
    
    itemprobs.group.ind <- rep(paste("Group", 1:G), times = sum(ncat[ which(model.indicator == 1) ]))
    itemprobs.var.ind <- rep(names(x$itemprob), times = G * ncat[ which(model.indicator == 1) ])
    itemprobs.cat.ind <- paste("Cat", rep(as.numeric(unlist(apply(t(ncat[ which(model.indicator == 1) ]), 2, function(x) 0:(x-1)))), each = G))
    
    x$classprob.initial <- xx$classprob
    x$itemprob.initial <- xx$itemprob
    
    x$samples <- list()
    
    if( MAP ) x$samples$logpost <- w$logpost else x$samples$loglik <- w$loglike
    x$samples$itemprob <- var.probs.l
    x$samples$classprob <- t( weights.mat )
    
    if(any(ncat > 2)){
      x$itemprob.tidy <- data.frame(itemprob = vec.itemprobs, group = itemprobs.group.ind, variable = itemprobs.var.ind, category = itemprobs.cat.ind)
    }else{ 
      # rearrange for backwards compatibility with plotting functions
      x$itemprob <- matrix(vec.itemprobs[itemprobs.cat.ind == "Cat 1"], nrow = G, ncol = sum(model.indicator), dimnames = list(paste("Group", 1:G), names(x$itemprob)))
      x$itemprob.sd <- matrix(vec.itemprobs.sd[itemprobs.cat.ind == "Cat 1"], nrow = G, ncol = sum(model.indicator), dimnames = list(paste("Group", 1:G), names(x$itemprob.sd)))
      M.in <- sum(model.indicator)
      arr <- array( dim=c(B,G,M.in) )
      for( m in 1:M.in ) 
      {
        arr[,,m] <- matrix( var.probs.l[[m]][,1], nrow=B , ncol=G , byrow = TRUE )
      }
      dimnames(arr)[[3]] <- names(v.probs$mean)
      x$samples$itemprob <- arr
    }
    
    npars <- G * sum( model.indicator * ( ncat - 1 ) ) + G - 1
    
    llike <-  blca.compute.log.post( X, ncat, x$classprob, x$itemprob, model.indicator, prior=NULL )
    x$Z <- predict.blca.pt.estimate( x )
    #llike <- lp$llike
    x$AIC <- 2*llike - 2 * npars
    x$BIC <- 2*llike - npars * log(N)
    
    #inputs
    x$B <- B
    if( MAP ) x$prior <- list( alpha=alpha, beta=beta, delta=delta )
    lp <- blca.compute.log.post( X, ncat, x$classprob, x$itemprob, model.indicator, prior=x$prior )
    if( MAP ) x$logpost <- lp else x$loglik <- lp
    x$counts.n <- counts.n
    x$ncat <- ncat
    x$model.indicator <- model.indicator
    x$MAP <- MAP 
    
    x$method <- "boot"
    
    x <- blca.return.order( x )
    
    if( any(ncat > 2) ) class(x) <- c("blca.boot", "blca.multicat", "blca") else class(x) <- c("blca.boot", "blca")
    
    return(x)	
    
  }
