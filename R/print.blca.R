print.blca <- function(x, ...){
  
  z<- list()
  
  btau<- x$classprob
  btheta<- x$itemprob
  btau.se<- x$classprob.sd
  btheta.se<- x$itemprob.sd
  
  names(btau) <- paste(rep("Group", length(btau)), 1:length(btau))
  if(!is.null(btau.se)) names(btau.se) <- names(btau)
  
  z$ItemProb<- btheta
  z$MembershipProb<- btau
  if(!is.null(btheta.se)) z$se$ItemProb<- btheta.se
  if(!is.null(btau.se)) z$se$MembershipProb<- btau.se
  
  k <- match( x$method, c("em", "gibbs", "boot", "vb", "collapsed") )
  str.par.est.vec <- c("MAP Estimates", "Posterior Mean Estimates", "MAP Estimates", rep("Posterior Mean Estimates",2))
  if( k %in% c(1,3) )
  {
    if( !x$MAP )
      str.par.est <- "Maximum Likelihood Estimates"
    else
      str.par.est <- str.par.est.vec[k]
  }else{
    str.par.est <- str.par.est.vec[k]
  }
  
  cat("\n",str.par.est,":\n",sep="")
  cat("\nItem Probabilities:\n", "\n",sep="")
  if( !any(x$ncat > 2) )
  {
    print(round(z$ItemProb,3))
  }else{
    nam <- names( z$ItemProb ) 
    for( l in 1:length(nam) )
    {
      cat("\n",nam[l],"\n")
      print( round(z$ItemProb[[l]], 3) )
    }
  }
  
  cat("\nClass Probabilities:\n", "\n",sep="")
  print(round(z$MembershipProb,3))
  
  if(!is.null(btheta.se)){	
    
    str.par.unc.vec <- rep( "Posterior Standard Deviation Estimates", 5)
    if( k %in% c(1,3) )
    {
      if(!x$MAP) 
        str.unc.est <- "Standard Error Estimates" 
      else
        str.unc.est <- str.par.unc.vec[k]
    }else{
      str.unc.est <- str.par.unc.vec[k]
    }
    
    cat("\n\n", str.unc.est ,":\n", sep="")
    cat("\nItem Probabilities:\n", "\n",sep="")
    if( !any(x$ncat > 2) )
    {
      print(round(z$se$ItemProb,3))
    }else{
      nam <- names( z$se$ItemProb ) 
      for( k in 1:length(nam) )
      {
        cat("\n",nam[k],"\n")
        print( round(z$se$ItemProb[[k]], 3) )
      }
    }
    
    cat("\nMembership Probabilities:\n", "\n")
    
    print(round(z$se$MembershipProb, 3))
  } else warning("Posterior standard deviations not returned.", call.=FALSE)
  
}
