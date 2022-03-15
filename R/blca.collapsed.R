blca.collapsed <- function( X, G, formula=NULL, ncat=NULL, alpha=1, beta=1, delta=1, 
                            start.vals=c("single"), counts.n=NULL, iter=5000, burn.in=500, thin=1, 
                            G.sel=TRUE, var.sel=FALSE, post.hoc.run=TRUE, control.post.hoc=list(iter=2000,burn.in=500,thin=1), 
                            var.prob.thresh=0.75, n.gibbs=nrow(X), only.gibbs=TRUE, G.max=30, 
                            G.prior=dpois(1:G.max, lambda=1)/sum(dpois(1:G.max, lambda=1)), # normalize
                            prob.inc=0.5, hprior.model=FALSE, relabel=TRUE, verbose=FALSE, verbose.update=1000 )
{
  
  args.passed <- as.list( environment() )
  #list of returns
  x <- list()
  x$call <- match.call()
  x$args <- args.passed
  
  # convert X into a numeric matrix for passage to C (index from 0) and check inputs
  D <- blca.prep.data( X, formula, counts.n, ncat )
  
  X <- D$X
  ncat <- D$ncat
  x$args$formula <- D$formula
  x$args$ncat <- ncat
  
  if( !is.null(D$missing.idx) )
  {
    warning("Missing values encountered in X: rows with NA have been removed. Imputation is available for method = 'gibbs'.", call.=FALSE)
    X <- na.omit(X)	  
  }
  
  N <- nrow(X)
  M <- ncol(X)
  
  if( G.sel == FALSE ) G.max <- G
  
  out.prior <- blca.check.prior.collapsed( alpha, beta, delta, G.sel, G, G.max, M, ncat, only.gibbs )
  prior.init.type <- out.prior$prior.init.type
  gamma <- out.prior$gamma
  delta <- out.prior$delta
  only.gibbs <- out.prior$only.gibbs
  
  if( G.sel & G.max < G ){
    stop("the maximum number of groups cannot be less than the initial number of groups: rerun with G.max (default 30) appropriately specified")
  }
  
  # thinning parameter
  if( thin > 1 || thin == 0 ) stop("argument thin gives a thinning rate and should be less than or equal to 1")
  
  Thinby <- floor(1/thin)
  
  stored <- floor( iter / Thinby )	
  
  memberships <- numeric(stored*N)
  variable.inclusion.indicator <- numeric(stored*M)
  num.groups <- numeric(stored)
  log.post <- numeric(stored)
  prior.include <- numeric(stored)
  
  #if( is.null(model.indicator) ) 
  model.indicator <- rep(1,M) 
  
  hparam <- c( delta[1], beta[1] )
  accrts <- numeric(6)
  
  w <- .C(	"BLCA_VS",																    as.integer(t(X)),
           as.integer(N),															  	as.integer(M),
           as.integer(ncat),															as.double(hparam),
           as.integer(prior.init.type),									  as.double(delta),
           as.double(gamma),
           as.integer(!G.sel),														as.integer(only.gibbs),
           as.integer(G),																  as.integer(G.max),
           as.integer(iter),															as.integer(burn.in),
           as.integer(Thinby),														as.integer(n.gibbs),
           memberships = as.integer(memberships),				  ngroups = as.integer(num.groups),
           as.double(G.prior),														as.integer(var.sel),
           variable.inclusion.indicator = as.integer(variable.inclusion.indicator),	as.double(prob.inc),
           log.posterior = as.double(log.post),					  as.integer(hprior.model),
           prior.include = as.double(prior.include),			as.integer(model.indicator),
           as.integer(verbose),													  as.integer(verbose.update),
           accrts=as.double(accrts),
           PACKAGE = "BayesLCA" )
  
  
  if( verbose ) cat("Finished sampling...\n")
  
  membership.mat <- matrix( w$memberships, nrow = stored, ncol=N, byrow=FALSE ) + 1
  vindicator.mat <- matrix( w$variable.inclusion.indicator, nrow=stored, ncol=M, byrow=FALSE )
  
  rl.flag <- ( G.sel == FALSE & G == 1 )
  
  if( relabel ) 
  {
    # relabel by the number of groups
    grp <- sort(unique(w$ngroups))
    idxs <- vector(length(grp),mode="list") 
    labels <- vector(length(grp),mode="list")
    varindicator <- vector(length(grp),mode="list")
    logpost <- vector(length(grp),mode="list")
    
    for( k in seq_along(grp) ) 
    {
      idxs[[k]] <- which( w$ngroups == grp[k] )
      idx.pivot <- idxs[[k]][ which.max( w$log.posterior[idxs[[k]]] ) ]
      mi.memb.mat <- membership.mat[idxs[[k]], , drop=FALSE]
      z.pivot <- membership.mat[ idx.pivot, ]
      s.z.pivot <- sort( table(z.pivot), decreasing=T )
      ma <- as.numeric( names(s.z.pivot) )
      z.pivot <- match( z.pivot, ma ) # label largest group 1, second largest group 2 etc.
      # sort by group size for nice presentation and plotting
      # call ecr function to get optimal permutations
      if( !rl.flag )
      {
        rl <- ecr( zpivot=z.pivot, mi.memb.mat, K=grp[k])
        ro.membership.mat <- array( dim=dim(mi.memb.mat))
        for( j in 1:nrow(mi.memb.mat) ) ro.membership.mat[j,] <- pmatch(  mi.memb.mat[j,], rl$permutation[j,], duplicates.ok=T )
        labels[[k]] <- ro.membership.mat
      }else{
        labels[[k]] <- mi.memb.mat
      }
      varindicator[[k]] <- vindicator.mat[idxs[[k]],]
      logpost[[k]] <- w$log.posterior[ idxs[[k]] ]
    }
    names(labels) <- names(varindicator) <- names(logpost) <- paste0("G = ",grp)
    
    # extract the most visited number of groups
    t <- table(w$ngroups)
    j <- which.max(t)
    Ghat <- as.numeric(names(t)[j])
    
    Z <- t( apply( labels[[j]], 2, tabulate, nbins=Ghat ) ) / nrow(labels[[j]])
  }else{
    Z <- NULL
    if( G == 1 & G.sel == FALSE ) Z <- matrix( rep(1,N), nrow=N, ncol=1 ) 
  }
  
  #x = list()
  
  #x$call <- match.call()
  
  x$classprob <- NULL
  x$classprob.sd <- NULL
  x$itemprob <- NULL
  x$itemprob.sd <- NULL
  
  x$Z <- Z
  
  x$samples <- list( logpost=logpost, G=w$ngroups, labels=labels, var.ind=varindicator )
  acc <- w$accrts
  x$accrates <- list()
  if( !only.gibbs ) x$accrates$labels <- list( m1 = acc[1] , m2 = acc[2] , m3 = acc[3] )
  if( G.sel ) x$accrates$G <- list( eject = acc[4], absorb = acc[5] )
  if( var.sel ) x$accrates$var.sel <- acc[6]
  
  if( var.sel )
  {
    x$samples$var.ind <- varindicator
    if(!is.null(colnames(X))) 
    {
      for(k in seq_along(grp) ) colnames(x$samples$var.ind[[k]]) <- colnames(X)
    }
  }
  
  if( hprior.model ) x$samples$prob.inc <- w$prior.include
  
  
  #inputs
  x$iter <- iter
  x$burn.in <- burn.in
  x$thin <- thin
  x$G.sel <- G.sel
  x$var.sel <- var.sel
  x$prob.inc <- prob.inc
  x$hprior.model <- hprior.model
  x$model.indicator <- model.indicator
  x$relabel <- relabel
  x$prior <- list( alpha=alpha, beta=beta, delta=delta, prob.inc=prob.inc, G.prior=G.prior, G.max=G.max )
  x$ncat <- ncat
  
  if( post.hoc.run ){
    # use the most visited 'G' model and run auxiliary run for plotting
    if( var.sel )
    {
      Mstar <- which( colMeans(x$samples$var.ind[[j]]) > var.prob.thresh )
    }else{
      Mstar <- which( model.indicator > 0 )
    }
    mod.star <- rep(0,M)
    mod.star[Mstar] <- 1
    x$model.indicator <- mod.star
    if( verbose ) cat("Running post hoc Gibbs sampling on most probable model...\n")
    ph.est <- blca.gibbs( x$args$X, G=Ghat, formula=formula, model.indicator=mod.star, ncat=ncat, 
                          iter=control.post.hoc$iter, burn.in=control.post.hoc$burn.in, thin=control.post.hoc$thin, verbose=FALSE )
    #cat("\nFinished post hoc Gibbs sampling...\n")
    x$classprob <- ph.est$classprob
    x$classprob.sd <- ph.est$classprob.sd
    x$itemprob <- ph.est$itemprob 
    x$itemprob.sd <- ph.est$itemprob.sd 
    x$itemprob.tidy <- ph.est$itemprob.tidy
    x$Z <- ph.est$Z
    x$samples$classprob <- ph.est$samples$classprob
    x$samples$itemprob <- ph.est$samples$itemprob
    if(verbose) cat("Posterior group membership probabilities computed from post hoc Gibbs sampling\n")
  } 
  
  x$method <- "collapsed"
  
  x <- blca.return.order( x )
  
  if( any(ncat > 2 ) ) class(x) <- c("blca.collapsed", "blca.multicat", "blca") else class(x)<-c("blca.collapsed", "blca")
  
  return(x)
  
}
