blca.gibbs <- function( X, G, formula=NULL,  ncat=NULL,  alpha=1, beta=1, delta=1, start.vals=c("prior","single","across"), 
                        counts.n=NULL, model.indicator=NULL, impute.missing=FALSE, iter=5000, burn.in=100, thin=1, accept=thin, 
                        relabel=TRUE, verbose=TRUE, verbose.update=1000 ) 
{
  # check if data is simulated 
  #if( class(X)[1] == "blca.rand" & !is.matrix(X) ) X <- X$X

  #list of returns
  args.passed <- as.list( environment() )
  x <- list()
  x$call <- match.call()
  x$args <- args.passed
  
  # convert X into a numeric matrix and check inputs
  D <- blca.prep.data( X, formula, counts.n, ncat )
  
  X <- D$X
  ncat <- D$ncat
  x$args$formula <- D$formula
  x$args$ncat <- ncat
  levnames <- D$levnames
  
  if( !is.null(D$missing.idx) & impute.missing )
  { 
    missing <- TRUE
    missing.idx <- D$missing.idx
  }else{ 
    missing <- FALSE 
  }
  
  if( !is.null(D$missing.idx) & !impute.missing )
  {
    warning("Missing values encountered in X and rows with NA have been removed: imputation is performed if impute.missing is TRUE", call.=FALSE )
    X <- na.omit(X)
  }
  
  N<-nrow(X) 
  M<-ncol(X)
  
  model.indicator <- rep(1,M)
  M.in <- sum( model.indicator ) 
  if( M.in == 0 ) stop("there are no variables included: check model.indicator")
  if( prod( ncat[model.indicator==1] ) <= (M.in+1)*G) warning( "Gibbs sampling can be run but check constraints on parameter estimability", call.=FALSE )
  
  out.prior <- blca.check.prior( alpha, beta, delta, G, M, ncat )
  prior.init.type <- out.prior$prior.init.type
  gamma <- out.prior$gamma
  delta <- out.prior$delta
  
  if( is.character(start.vals[1]) ) 
  {
    if( start.vals[1] == "single" ) init <- 0 
    if( start.vals[1] == "across" ) init <- 1
    if( start.vals[1] == "prior" ) init <- 2
  }else{
    init <- 0
  }
  
  # thinning parameter
  if( thin > 1 || thin == 0 ) stop("argument thin gives a thinning rate and should be less than or equal to 1")
  # thinning parameter
  Thinby <- floor(1/thin)
  
  #if((iter-burn.in)%%thin != 0)
  #	stop("\t Please thin by an amount that divides perfectly into iter - burn.in. ")
  
  stored <- floor( iter / Thinby )
  
  memberships <- numeric(stored*N)
  weights <- numeric(stored*G)
  variable.probs <- numeric(stored*G*sum(ncat))
  log.post <- numeric( stored )
  log.like <- numeric( stored )
  
  # set up storage for imputed missing data if present
  
  if( missing & impute.missing )
  {
    #warning("encountered missing values in X: these will be imputed using sampling")
    # add an na behaviour option in time
    n.missing <- nrow(missing.idx)
    missing.values <- numeric( n.missing * stored )
    # the position of the missing values in C ordering
    #  pass the row and the corresponding column, as easier to find... 
    sample.missing.data <- TRUE
    # replace the missing values with random generated
    for( k in 1:n.missing ) X[ missing.idx[k,1], missing.idx[k,2] ] <- sample(0:(ncat[missing.idx[k,2]]-1), size=1)
    position.missing <- as.vector(t(missing.idx)) - 1
  }else{
    n.missing <- 0
    missing.values <- NULL
    position.missing <- NULL
    sample.missing.data <- FALSE
  }
  
  hparam <- c( delta[1], beta[1] )
  
  # initialization is done randomly to groups
  x.vec <- as.vector(t(X))
  
  w <- .C(	"BLCA_GIBBS_SAMPLER",											as.integer(x.vec),
           as.integer(N),																as.integer(M),
           as.integer(ncat),															as.double(hparam),
           as.integer(prior.init.type),
           as.double(delta),															as.double(gamma),
           as.integer(init),
           as.integer(G),															 	as.integer(iter),											
           as.integer(burn.in),
           as.integer(Thinby),														memberships = as.integer(memberships),					
           weights = as.double(weights),									variable.probs = as.double(variable.probs),			
           as.integer(model.indicator), 									log.posterior = as.double(log.post), 
           log.like = as.double(log.like),
           as.integer(sample.missing.data),                as.integer(n.missing),
           missing.values = as.integer(missing.values),    as.integer(position.missing),
           as.integer( verbose ),													as.integer( verbose.update ),			
           PACKAGE = "BayesLCA" )
  
  if( verbose ) cat("\nFinished sampling...\n")
  
  membership.mat <- matrix( w$memberships, nrow = stored, ncol=N, byrow=TRUE ) + 1
  
  #split these into a list as it will be easier to rearrange from label processing
  var.probs.l <- vector( length=sum(model.indicator), mode="list")
  
  l <- 1
  for( j in 1:M )
  {
    if(j == 1){
      gap <- 0
    }else{
      gap <- stored * G * sum( ncat[1:(j-1)] )
    }
    #variable probabilities stacked by group and then iteration
    if( model.indicator[j] )
    {
      var.probs.l[[l]] <- matrix( w$variable.probs[(gap+1):(gap + stored*G*ncat[j])] , nrow = stored * G, ncol=ncat[j], byrow=TRUE )
      l <- l+1
    }
  }	
  
  weights.mat <- matrix( w$weights , nrow=stored, ncol=G, byrow=TRUE )
  
  if( relabel ) 
  {
    if( verbose ) cat("\nPost processing samples...\n")
    
    # get the pivoting label vector
    idx.pivot <- which.max( w$log.like )
    z.pivot <- membership.mat[ idx.pivot, ]
    # sort by group size for nice presentation and plotting
    #s.z.pivot <- sort( table(z.pivot), decreasing=T )
    #ma <- as.numeric( names(s.z.pivot) )
    #z.pivot <- match( z.pivot, ma )
    
    # call ecr function to get optimal permutations
    rl <- ecr( zpivot = z.pivot, membership.mat, K=G)
    
    #re-order weights and probabilities according to ecr permutations
    perm <- rl$permutations
    for( k in 1:stored )
    {
      # this does not appear to be working so well
      weights.mat[ k, ] <- weights.mat[ k, perm[k,] ]
      l <- 1
      for( j in 1:M  )
      {
        it <- (k-1)*G
        if( model.indicator[j] )
        {
          var.probs.l[[l]][ it + 1:G , ] = var.probs.l[[l]][ it + perm[ k,] , ]
          l <- l+1
        }
      }
    }
    
    # re-order membership matrix acccording to permutations
    ro.membership.mat <- array( dim=dim(membership.mat) )
    for( k in 1:stored )
    {
      ro.membership.mat[k,] <- pmatch( membership.mat[k,], perm[k,], duplicates.ok=T )
    }
    
    membership.mat <- ro.membership.mat
    
    Z <- t( apply( membership.mat, 2, tabulate, nbins=G ) ) / stored  ## re-do Z using predict?
    
  }else{
    Z <- NULL
  }
  
  x$classprob <- apply( weights.mat, 2, mean )
  x$classprob.sd <- apply( weights.mat, 2, sd )
  o <- order(x$classprob, decreasing=TRUE)
  x$classprob <- x$classprob[o]
  x$classprob.sd <- x$classprob.sd[o]
  Z <- Z[,o]
  
  for( k in 1:stored )
  {
    membership.mat[k,] <- pmatch( membership.mat[k,], o, duplicates.ok=T )
  }
  
  l.len <- sum(model.indicator)
  v.probs <- list( mean=vector(l.len,mode="list"), sd=vector(l.len,mode="list") )
  tt <- seq( 0, G*(stored-1), by=G ) # this should be G*stored - 1?
  l <- 1
  for( j in 1:M )
  {
    if( model.indicator[j] )
    {
      v.probs$mean[[l]]  <- matrix( nrow = G, ncol = ncat[j] ) 
      v.probs$sd[[l]] <- matrix( nrow = G, ncol = ncat[j] )
      for( g in 1:G )
      {
        b <- var.probs.l[[l]][ tt+g , ]
        v.probs$mean[[l]][g,] <- apply( b, 2, mean )
        v.probs$sd[[l]][g,] <- apply( b, 2, sd )
      }
      v.probs$mean[[l]] <- v.probs$mean[[l]][o,]
      v.probs$sd[[l]] <- v.probs$sd[[l]][o,]
      rownames( v.probs$mean[[l]] ) <- paste( "Group", 1:G )
      colnames( v.probs$mean[[l]] ) <- levnames[[j]] 
      rownames( v.probs$sd[[l]] ) <- paste( "Group", 1:G )
      colnames( v.probs$sd[[l]] ) <- levnames[[j]] 
      l <- l+1
    }
  }
  
  # finally reorder the samples for correct prediction of new data 
  weights.mat <- weights.mat[ , o]
  for( k in 1:stored )
  {
    l <- 1
    for( j in 1:M  )
    {
      it <- (k-1)*G
      if( model.indicator[j] )
      {
        var.probs.l[[l]][ it + 1:G , ] = var.probs.l[[l]][ it + o , ] # final reorder as label swapping done
        l <- l+1
      }
    }
  }
  
  
  if( is.null(colnames(X)) )
  {
    names( v.probs$mean ) <- names( v.probs$sd ) <- paste("Variable", which( model.indicator==1 ) )
  }else{
    names( v.probs$mean ) <- names( v.probs$sd ) <- colnames(X)[  which( model.indicator == 1 ) ]
  }
  
  if( missing )
  {
    missing.samp <- matrix( w$missing, nrow=stored, byrow=TRUE ) + 1
    missing.samp.lev <- array( NA, dim=dim(missing.samp) )
    # match levels to sampled cats for presentation
    vars.missing <- missing.idx[,2] 
    for( k in seq_along(vars.missing) ) # k indexes cols of missing.samp
    {
      j <- vars.missing[k]
      missing.samp.lev[,k] <- levnames[[j]][ missing.samp[,k] ] 
    }
    missing.samp <- missing.samp.lev
  }
  
  # compile the list to return
  #x$classprob <- apply( weights.mat, 2, mean )
  
  x$itemprob <- v.probs$mean
  
  #x$classprob.sd <- apply( weights.mat, 2, sd )
  x$itemprob.sd <- v.probs$sd
  
  x$Z <- Z
  x$logpost <- max(w$log.posterior)
  
  # look after itemprob.tidy here
  vec.itemprobs <- unlist(x$itemprob)
  vec.itemprobs.sd <- unlist(x$itemprob.sd)
  
  itemprobs.group.ind <- rep(paste("Group", 1:G), times = sum(ncat[ which(model.indicator == 1) ]))
  itemprobs.var.ind <- rep(names(x$itemprob), times = G * ncat[ which(model.indicator == 1) ])
  # need to rewrite this with the new names
  # extract elements from levnames for which model.indicator == 1
  newlevnames <- vector( length=sum(model.indicator), mode="list")
  c <- 1
  for( k in which( model.indicator == 1 ) )
  {
    newlevnames[[c]] <- levnames[[k]]
    c <- c + 1
  }
  z <- unlist( newlevnames )
  itemprobs.cat.ind <- rep.int( z, times=rep(G,length(z)) )
  
  x$samples <- list()
  
  x$samples$logpost <- w$log.posterior
  x$samples$loglik <- w$log.like
  x$samples$Giter <- rep( G, stored )
  x$samples$itemprob <- var.probs.l
  x$samples$classprob <- weights.mat 
  x$samples$labels <- membership.mat
  
  if( missing ) x$samples$missing <- missing.samp
  
  if(any(ncat > 2)){
    x$itemprob.tidy <- data.frame(itemprob = vec.itemprobs, group = itemprobs.group.ind, variable = itemprobs.var.ind, category = itemprobs.cat.ind)
  }else{ 
    # rearrange for backwards compatibility with plotting functions
    
    # need to be careful here as the levels have been rearranged - take second entry of newlevnames
    idx <- NULL
    for( k in seq_along(newlevnames) ) idx <- union( idx, which( itemprobs.cat.ind == newlevnames[[k]][2] ) )
    x$itemprob <- matrix(vec.itemprobs[ idx ], nrow = G, ncol = sum(model.indicator), dimnames = list(paste("Group", 1:G), names(x$itemprob)))
    x$itemprob.sd <- matrix(vec.itemprobs.sd[ idx ], nrow = G, ncol = sum(model.indicator), dimnames = list(paste("Group", 1:G), colnames(x$itemprob)))
    
    M.in <- sum(model.indicator)
    arr <- array( dim=c(stored,G,M.in) )
    for( m in 1:M.in ) 
    {
      arr[,,m] <- matrix( var.probs.l[[m]][,2], nrow=stored , ncol=G , byrow = TRUE )
    }
    dimnames(arr)[[3]] <- names(v.probs$mean)
    #x$samples$itemprob <- arr
  }
  
  # for information criteria take the mean of log-likelihoods 
  Dbar<- mean(w$log.like)
  S2<- var(w$log.like)
  
  # now get the log-likelihood at the MCMC Bayes estimate of parameters
  log.like.val <- 0.0
  weights.be <- x$classprob
  variable.probs.be <- unlist( lapply(x$itemprob, t) )
  
  w.ll <- .C(	"BLCA_LOG_LIKE",											        as.integer(x.vec),
              as.integer(N),																as.integer(M),
              as.integer(ncat),															as.double(hparam),
              as.integer(G),															 				
              weights = as.double(weights.be),							variable.probs = as.double(variable.probs.be),			
              as.integer(model.indicator), 									value = as.double(log.like.val),
              PACKAGE = "BayesLCA" )
  
  loglik.be <- w.ll$value
  
  # DIC is computed correctly here?  
  #   x$samples$logpost has to be replaced by the log posterior evaluated i.e. x$logpost
  x$DIC<- 2*(2*Dbar - loglik.be)
  # need to check definition of BICM (number of samples)
  x$BICM<- 2*(Dbar - S2*(log(N)-1))
  x$AICM<- 2*(Dbar - S2)
  
  #inputs
  x$G <- G
  x$iter <- iter
  x$burn.in <- burn.in
  x$thin <- thin
  x$relabel <- relabel
  x$prior <- list( alpha=alpha, beta=beta, delta=delta )
  x$ncat <- ncat
  x$model.indicator <- model.indicator
  
  x$method <- "gibbs"
  
  x <- blca.return.order( x )
  
  if( any(ncat > 2) ) class(x) <- c("blca.gibbs", "blca.multicat", "blca") else class(x) <- c("blca.gibbs", "blca")
  
  return(x)
}

