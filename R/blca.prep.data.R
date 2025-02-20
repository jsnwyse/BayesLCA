blca.prep.data <- function( X, formula=NULL, counts.n, ncat )
{
  # check if data is simulated 
  if( class(X)[1] == "blca.rand" & !is.matrix(X) ) X <- X$X
  # check if data is "data.blca" type
  if( class(X)[1] == "data.blca" || !is.null(counts.n) )
  {
    if( class(X)[1] == "data.blca" ) 
    {
      z <- X$counts.n
      Y <- X$data 
    }else{ 
      z <- counts.n
      Y <- as.matrix(X)
    }
    U <- matrix( rep( Y[1,], z[1] ) , nrow=z[1], byrow=TRUE )
    for( k in 2:length(z) ) U <- rbind( U, matrix( rep( Y[k,], z[k] ) , nrow=z[k], byrow=TRUE ) )
    X <- as.matrix(U)
  }
  if( class(X)[1] == "blca.rand" ) X <- unclass(X)
  
  if( is.matrix(X) ) X <- as.data.frame( X )
  
  # create formula for prediction methods
  if( is.null(formula) ) formula <- as.formula(paste0( "~", paste( names(X), collapse = "+")))
  X <- model.frame( formula, data=X, na.action = na.pass )
  
  # check for missing values
  missing <- any(is.na(X))
  if( missing ) idxs <- which( is.na(X), arr.ind=T ) else idxs <- NULL
  
  # check here if X is a data frame
  N <- nrow(X)
  M <- ncol(X)
  
  Y <- matrix( nrow=N, ncol=M )
  colnames(Y) <- colnames(X) 
  
  levnames <- vector( M, mode="list" )
  
  # convert everything to factors and match Y vals
  for( k in 1:M )
  {
    X[,k] <- factor( X[,k] )
    levnames[[k]] <- levels( X[,k] )
    Y[,k] <- match( X[,k], table = levnames[[k]] )
  }
  
  #check that matrix is binary and/or ncat is passed
  if( is.null(ncat) )
  {
    # determine ncat from the levnames vector
    ncat <- unlist( lapply( levnames, length ) )
  }else{
    if(length(ncat)!= M) stop("the length of ncat must be the same as the number of variables")
    if( any( ncat != unlist( lapply(levnames,length) ) ) ) stop("if passing ncat, the number of levels within variables must correspond to ncat entries")
  }
  
  # subtract 1 from Y to offset at zero
  return( list( X=Y-1, ncat=ncat, levnames = levnames, missing.idx = idxs, formula = formula ) )
}


