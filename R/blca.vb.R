blca.vb <- function( X, G, formula=NULL,  ncat=NULL, alpha=1, beta=1, delta=1, start.vals = c("single", "across"), 
            counts.n=NULL, iter=5000, restarts=5, verbose=TRUE, conv=1e-6, small=1e-10 )
  {
    
    #list of returns
    args.passed <- as.list( environment() )
    x <- list()
    x$call <- match.call()
    x$args <- args.passed
    
    # convert X into a numeric matrix and check inputs
    D <- blca.prep.data( X, formula, counts.n, ncat )
    
    X <- D$X
    x$args$formula <- D$formula
    ncat <- D$ncat
    x$args$ncat <- ncat
    levnames <- D$levnames
    x$G <- G
    
    if( !is.null(D$missing.idx) )
    {
      warning("Missing values encountered in X: rows with NA have been removed", call.=FALSE)
      X <- na.omit(X)	  
    }
    
    N<-nrow(X) 
    M<-ncol(X)
    
    model.indicator <- rep(1,M)
    M.in <- sum( model.indicator ) 
    if( M.in == 0 ) stop("there are no variables included: check model.indicator")

    out.prior <- blca.check.prior( alpha, beta, delta, G, M, ncat )
    prior.init.type <- out.prior$prior.init.type
    gamma <- out.prior$gamma
    delta <- out.prior$delta
        
    if(is.numeric(restarts)){
      if(length(restarts)>1){
        restarts<- restarts[1]
        warning("Restarts improperly specified: first value will be used, other values will be ignored", call.=FALSE)
      }# else {stop("restarts improperly specified. Must be an integer of length 1.")}
    } else {stop("restarts improperly specified: must be integer valued")}
    
    #storage 
    weights <- pars.weights <- numeric(G)
    sd.weights <- numeric(G)
    variable.probs <- pars.variable.probs <- numeric(G*sum(ncat))
    sd.variable.probs <- numeric(G*sum(ncat))
    group.probs <- numeric(N*G)
    lb <- numeric( iter )	
    iters <- 0
    converged <- 0
    logpost <- 0
    
    hparam <- c( delta[1], beta[1] )
    
    lb.max <- -Inf
    w.max <- list()
    
    init.vals <- numeric(N*G)
    if( is.character(start.vals[1]) ) 
    {
      if( start.vals[1] == "single" ) init <- 0 else init <- 1
    }else{
      init <- 2
      if(is.matrix(start.vals) & all(dim(as.matrix(start.vals)) == c(N,G))){
        init.vals<- as.vector(start.vals) 
      }else{
        if(is.numeric(start.vals) & length(as.numeric(start.vals))==N)
        { 
          init.vals<- as.vector(unMAP(start.vals))
        }else{ 
          stop("start.vals improperly specified: see help")
        }
      }
    }
    
    cnv.warn <- FALSE
    
    for( r in 1:restarts )
    {	
      # call the VB algorithm
      
      w <- .C( 		"BLCA_VB_FIT", 										as.integer(t(X)), 			
                 as.integer(N),											as.integer(M), 
                 as.integer(ncat),	 									as.double(hparam),
                 as.integer(prior.init.type),						as.double(delta),
                 as.double(gamma),
                 as.integer(G),											
                 as.integer(init), 									as.double(init.vals),
                 as.integer(iter),
                 iters = as.integer(iter),							group.probs = as.double(group.probs),
                 weights = as.double(weights),						se.weights = as.double( sd.weights ),
                 pars.weights = as.double(pars.weights),
                 variable.probs = as.double(variable.probs),	se.variable.probs = as.double( sd.variable.probs ),
                 pars.variable.probs = as.double(pars.variable.probs),
                 as.integer(model.indicator),						lb = as.double(lb),
                 as.double(conv),										converged = as.integer(converged),
                 logpost = as.double(logpost),
                 PACKAGE = "BayesLCA" )
      
      #need to extract the optimizer parameters for vb here...
      
      lb.this <- w$lb[ w$iters ]
      
      #store the run that gives the highest value  of the log posterior
      #	for the runs that have converged
      new.max <- FALSE
      if( lb.this > lb.max )
      {
        w.max <- w 
        lb.max <- lb.this
        new.max <- TRUE
      }
      
      if( verbose ) 
      {
        if( new.max && r>1 )
        {
          cat( "\nNew maximum found... Restart number ",r,", lower bound = ", round(lb.this,2),"...", sep = "" )
        }else{ 
          cat( "\nRestart number ",r,", lower bound = ", round(lb.this,2),"...", sep = "" )
        }
      }
      
      if( as.logical(w$converged) == FALSE ) cnv.warn <- TRUE 
      
    }	
    if( cnv.warn ) warning("Some restarts failed to converge: rerun with a higher iter value or less stringent tolerance", call.=FALSE)
    if( verbose ) cat("\n")
    
    w <- w.max
    
    # an issue here  with  ordering the itemprobs
    
    var.probs.l <- list()
    se.var.probs.l <- list()
    par.var.probs.l <- list()
    
    x$classprob <- w$weights
    o <- order( x$classprob, decreasing = TRUE )
    x$classprob <- x$classprob[o]
    
    x$classprob.sd <- w$se.weights[o]
    
    x$Z <- matrix( w.max$group.probs, nrow=N, ncol=G, byrow=TRUE )
    x$Z <- x$Z[,o]
    colnames(x$Z) <- paste( "Group", 1:G )
    
    l <- 1
    for( j in 1:M )
    {
      if(j == 1){
        gap <- 0
      }else{
        gap <- G * sum( ncat[1:(j-1)] )
      }
      #variable probabilities stacked by group 
      if( model.indicator[j] )
      {
        var.probs.l[[l]] <- matrix( w$variable.probs[(gap+1):(gap + G*ncat[j])] , nrow =  G, ncol=ncat[j], byrow=TRUE )
        var.probs.l[[l]] <- var.probs.l[[l]][o,]
        se.var.probs.l[[l]] <- matrix( w$se.variable.probs[(gap+1):(gap + G*ncat[j])] , nrow =  G, ncol=ncat[j], byrow=TRUE )
        se.var.probs.l[[l]] <- se.var.probs.l[[l]][o,]
        par.var.probs.l[[l]] <- matrix( w$pars.variable.probs[(gap+1):(gap + G*ncat[j])], nrow=G, ncol=ncat[j], byrow=TRUE )
        par.var.probs.l[[l]] <- par.var.probs.l[[l]][o,]
        rownames( var.probs.l[[l]] ) <- rownames( se.var.probs.l[[l]] ) <- rownames(par.var.probs.l[[l]]) <-  paste( "Group", 1:G )
        colnames( var.probs.l[[l]] ) <- colnames( se.var.probs.l[[l]] ) <- colnames(par.var.probs.l[[l]]) <-  levnames[[j]] #paste("Cat",0:(ncat[j]-1) )
        l <- l+1
      }
    }
    
    if( is.null(colnames(X)) )
    {
      names( var.probs.l ) <- names( se.var.probs.l ) <- paste("Variable", which( model.indicator==1 ) )
    }else{
      names( var.probs.l ) <- names( se.var.probs.l ) <- colnames(X)[  which( model.indicator == 1 ) ]
    }
    
    x$itemprob <- var.probs.l	
    #x$classprob <- w$weights
    
    x$itemprob.sd <- se.var.probs.l	
    #x$classprob.sd <- w$se.weights
    
    vec.itemprobs <- unlist( x$itemprob )
    vec.itemprobs.se <- unlist( x$itemprob.sd )
    
    x$parameters <- list( classprob=w$pars.weights[o], itemprob=par.var.probs.l )
    #x$parameters$classprob <- w$pars.weights
    #x$parameters$itemprob <- par.var.probs.l
    
    itemprobs.group.ind <- rep(paste("Group", 1:G), times = sum(ncat[ which(model.indicator == 1) ]))
    itemprobs.var.ind <- rep(names(x$itemprob), times = G * ncat[ which(model.indicator == 1) ])
    itemprobs.cat.ind <- paste("Cat", rep(as.numeric(unlist(apply(t(ncat[ which(model.indicator == 1) ]), 2, function(x) 0:(x-1)))), each = G))
    
    if( any(ncat > 2) ){
      x$itemprob.tidy <- data.frame(itemprob = vec.itemprobs, group = itemprobs.group.ind, variable = itemprobs.var.ind, category = itemprobs.cat.ind)
    }else{ 
      x$itemprob <- matrix(vec.itemprobs[itemprobs.cat.ind == "Cat 1"], nrow = G, ncol = sum(model.indicator), dimnames = list(paste("Group", 1:G), names(x$itemprob)))
      x$itemprob.sd <- matrix( vec.itemprobs.se[ itemprobs.cat.ind == "Cat 1" ], nrow=G, ncol=sum(model.indicator), dimnames = list(paste("Group", 1:G), names(x$itemprob)) )
      
      M.in <- sum(model.indicator)
      arr <- array(0,dim=c(G, M.in ,2) )
      for( j in 1:M.in ) 
      {
        arr[ , j , ] <- x$parameters$itemprob[[j]]
      }
      dimnames(arr)[[2]] <- names(var.probs.l)
      x$parameters$itemprob <- arr
    }
    
    #x$Z <- matrix( w.max$group.probs, nrow=N, ncol=G, byrow=TRUE )
    #colnames(x$Z) <- paste( "Group", 1:G )
    
    x$LB <- lb.max
    x$lbstore <- w$lb[ 1:w$iters ]
    
    
    x$converged <- as.logical(w$converged)
    x$iter <- w$iters
    x$eps <- w$lb[ w$iters ] - w$lb[ w$iters-1 ]
    
    x$prior<- list( alpha=alpha, beta=beta, delta=delta )
    x$model.indicator <- model.indicator
    x$ncat <-  ncat
    
    x$method <- "vb"
    
    x <- blca.return.order( x )
    
    if( any(ncat>2) ) class(x)<-c("blca.vb", "blca.multicat", "blca") else class(x) <- c("blca.vb", "blca" )
    
    return(x)
  }
