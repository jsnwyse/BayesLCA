plot.blca.multicat <- function(x, main="", variables = NULL, col1=NULL, ...)
{
  #logistore<- devAskNewPage()
  which.pl <-1L
  
  if( !is.null(col1) & any(x$args$ncat != x$args$ncat[1]) ) stop("can only use argument col1 when number of categories equal over all variables")
  
  if( !is.null(col1) & length(col1) != x$args$ncat[1] ) stop("col1 must have length corresponding to the number of categories")
  
  show <- rep(FALSE, 11)
  show[which.pl] <-TRUE
  
  if(show[1])
  {
    # may need to split the plotting window
    old.mfrow <- par()$mfrow
    
    tau <- x$classprob
    #itempr <- x$itemprob
    
    G <- length(tau)
    names(tau)<- paste("Group", 1:G)
    sepG <- 0.075/G
    
    # if subset of variables not specified, plot all
    if( is.null(variables) ) variables <- 1:length(x$itemprob)
    M <- length(variables)
    if( M > length(x$itemprob) ) stop("more variables specified than available")
    
    # do a maximum of four variable plots in a panel-- determine layout
    r <- min(M,4)
    if( r > 2 ) new.mfrow <- c(2,2) else new.mfrow <- c(1,2)
    # only reset if r > 1
    if( r > 1 ) oldpar <- par( mfrow=new.mfrow )
    
    if(M > 4) devAskNewPage(TRUE)
    
    for(v in variables)
    {
      # pick out variable probabilities
      theta <- x$itemprob[[v]]
      #determine splits
      ncat <- ncol(theta)
      sepcat<- 0.045/ncat
      totcat <- (ncat + 2) * sepcat + 1
      totgrp <- (G+2)*sepG + 1
      # sep up plotting window
      plot.new()
      plot.window( xlim=c(-sepG,totgrp), ylim=c(-sepcat,totcat), main = names(theta[1]))
      xleft <- rep.int( c(0,cumsum(tau)[1:(G-1)]) + 0:(G-1) * sepG, times=rep(ncat,G) )
      xright <- rep( cumsum(tau) + 0:(G-1) * sepG, times = rep(ncat,G) )
      xmid <- 0.5 * ( xleft + xright )[1+0:(G-1)*(ncat)]
      p <- theta[1,]
      p0 <- apply( theta, 1, cumsum )
      q0 <- rbind( rep(0,G), p0[1:(ncat-1),])
      ybottom <- as.vector(q0) + rep( 0:(ncat-1) * sepcat, G )  
      ytop <- as.vector(p0) + rep( 0:(ncat-1)*sepcat, G) 
      ymid <- 0.5 * ( ybottom + ytop )[1:ncat]
      if( is.null(col1) ) col2 <- rainbow(ncat) else col2 <- col1
      rect( xleft, ybottom, xright, ytop, density=23, col=col2 )  
      axis( 1, at = xmid, labels = rownames(theta), lwd=-1, cex.axis=1, las=2)
      axis( 2, at = ymid, labels = colnames(theta),lwd=-1, cex.axis=1, las=2)
      mstr <- names(x$itemprob)[v]
      if( is.null(mstr) ) mstr <- paste("Variable",v)
      mtext( mstr, 3, 0.25, cex=1)
      
      #if( r > 2) devAskNewPage()
    }
    
    if( r > 1 ) par(mfrow=c(1,1))
  }
  #Show2 Multicat Parameter Mosaicplot
  #old_mfrow <- par()$mfrow
  #show[1]<- FALSE
  #which <- which(show)
  if(M > 4) devAskNewPage(FALSE)
  #plot.blca(x, which, main = main, col1 = col1, ...)
}

