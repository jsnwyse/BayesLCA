undo.label.switching <- function( Z, Gsamp = NULL, logpost )
#undo label switching using label.switching package
{

	# Z is a matrix of labels with row i indexing classifications
  	# to groups from 1 : ngroups[i] or 1:ngroups if ngroups is an integer 
  	
  	ngroups <- Gsamp
  	
  	if( is.null(ngroups) ) stop("\t argument ngroups must be specified as either a fixed number of groups or a vector of the number of groups corresponding to each row of Z")
    
    if( length(ngroups) == 1 ) ngroups <- rep(ngroups, nrow(Z))
    
    #different no.'s groups
    
    Zrelab <- matrix( nrow = nrow(Z), ncol = ncol(Z) )
    
    Perm <- matrix( nrow = nrow(Z), ncol = max(ngroups) )
    
    nobs <- ncol(Z)
    
    G <- unique(ngroups)
    
    ret <- list()
    ret$groups <- ngroups
    ret$ncomponents <- numeric(length(G))
    ret$item.tags <- list()
    
    groups.not.done <- NULL
    
    j <- 1
    
    for(k in G)
    {
      
      idx.k <- which(ngroups == k)
      labels.k <- Z[idx.k,]
      logpost.k <- logpost[ idx.k ]
      
      idx.lab.k <- which.max(logpost.k)
      z.pivot <- labels.k[ idx.lab.k, ]
      
      nsamp.k <- length(idx.k)
      ngrp.k <- k
      
      permutation <- numeric(nsamp.k*k)
      
      if( length(idx.k) == 1 ) groups.not.done <- c( groups.not.done, k )

      if((k!=1) && (length(idx.k) > 1))
      {
      
        #nonempty.k <- apply( labels.k, 1, function(x) length(unique(x)) )
        #t.k <- sort( nonempty.k, index.return=TRUE )
        
        #labels.arranged.k <- labels.k[ t.k$ix, ]
        
        item.tags.k <- idx.k #[ t.k$ix ] #this is the actual index of each row in the original matrix Z
        
        #labels.out.k <- numeric(nsamp.k*nobs)
        
        #w <- .C(	"BLCA_RELABEL",								as.integer(nobs),
        #				as.integer(nsamp.k),					as.integer(ngrp.k),
        #				as.integer(labels.arranged.k),	x=as.integer(labels.out.k),
        #			   xx = as.integer(permutation),		PACKAGE = "BayesLCA" )
        
        w <- ecr( zpivot=z.pivot, labels.k, K=k)
        
        ret$ncomponents[j] <- k
        #ret$memberships[[j]] <- matrix(w$x,nrow = nsamp.k,ncol=nobs,byrow=FALSE)
        ret$permutation[[j]] <- w$permutations #matrix(w$xx,nrow=nsamp.k,ncol=k,byrow=FALSE)
        
        labels.perm <- labels.k
        # have to make a permuted membership matrix
        for( i in 1:nsamp.k )
        {
          labels.perm[i,] <- pmatch( labels.k[i,], w$permutations[i,], duplicates.ok=T )
        }
        ret$memberships[[j]] <- labels.perm
        
        #compute membership probabilities for each data pt
        probs <- matrix( nrow = nobs, ncol=k )
        
        for( id in 1:nobs ) probs[ id, ] <- tabulate( ret$memberships[[j]][,id], nbins=k )
        
        # for(id in 1:nobs){
        # 	for(c in 1:k){
        # 		probs[id,c] <- length(which(ret$memberships[[j]][,id] == c))
        # 	}
        # }
        # 
        
        probs <- probs/nsamp.k
        
        ret$membership.probabilities[[j]] <- probs
 
        #for variable indicator
		    ret$item.tags[[j]] <- idx.k  
		
	    	#store in the new Z matrix Zrelab
		
	    	Zrelab[ item.tags.k, ] <- ret$memberships[[j]]
		
	   	  Perm[ item.tags.k, 1:k ] <- ret$permutation[[j]]
        
      }else{
        
        ret$ncomponents[j] <- k
        ret$memberships[[j]] <- labels.k
        ret$membership.probabilities[[j]] <- matrix( 0 , nrow=nobs, ncol=k )
        for( id in 1:nobs ) ret$membership.probabilities[[j]][ id, labels.k[id] ] <- 1
        idx.k <- which(ngroups == k)
        ret$item.tags[[j]] <- idx.k #only in this case
        
      }
      
      j = j+1
      
    }
    
    
    ord <- sort( ret$ncomponents, index.return=TRUE )$ix
    
    membprob <- list()
    
   for( k in 1:length(ord) ) membprob[[k]] <- ret$membership.probabilities[[ ord[k] ]]
    
    names( membprob ) <- paste0("G=", sort(ret$ncomponents) ) 
    
    x <- list()
    
    x$call <- match.call()
    
    x$relab <- Zrelab
    
    x$numcomponents <- sort( ret$ncomponents )
    
    x$components <- sort( ret$ncomponents )
    
    x$label.probs <- membprob
    
    x$permutation <- Perm
    
    return(x)
    
  
  
}

