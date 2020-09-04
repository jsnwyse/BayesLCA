# plot.blca.boot <-
# function(x, which=1L, main="", ...){
#   show<- rep(FALSE, 5)
#   show[which]<- TRUE
#   if(show[5]){ 
#   	warning("no diagnostic plot available for boot method")
#   	show[5]<- FALSE
#   	}
#   which<- c(1:4)[show]
#   NextMethod("plot")
# }
