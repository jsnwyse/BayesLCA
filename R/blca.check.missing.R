blca.check.missing <- function( X )
{
  # scan the data frame/matrix for missing data
  missing <- any(is.na(X))
  
  if( missing )
  {
    idxs <- which( is.na(X), arr.ind=T )
  }else{
    idxs <- NULL
  }
  
  return(list( missing=missing, idxs=idxs ))
}
