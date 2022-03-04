blca.predict.check.data <- function( olddata, newdata=NULL )
{
  if( is.null(newdata) ) newdata <- olddata
  
  if( attributes(newdata)$class != attributes(olddata)$class ) 
    stop("newdata is not of class ",attributes(olddata)$class)
  
  if( !identical( newdata$names, olddata$names ) ){}
  # need to resolve naming here... 
}