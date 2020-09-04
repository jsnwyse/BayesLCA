plot.blca.em <-
function(x, which=1L, main="", ...){
  show<- rep(FALSE, 5)
  show[which]<- TRUE
  if(any(show[3:4])){ 
  	warning("density plots not supported for EM fits")
  	show[3:4]<- FALSE
  	 }
  which<- c(1:4,7)[show]
  NextMethod("plot")
}
