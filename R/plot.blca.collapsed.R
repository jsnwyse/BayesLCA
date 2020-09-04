# plot.blca.collapsed <-
#   function(x, which=1L, main="", ...){
#     #class(x)<- "blca"
#     #print("NextMethodUsed")
#     show<- rep(FALSE, 5)
#     show[which]<- TRUE
#     if(show[1] & is.null(x$itemprob)){
#       warning("item probabilities not calculated")
#       show[1]<- FALSE 
#     }
#     if(show[1] & ( x$G.sel | x$var.sel)) warning("optimal model from collapsed sampling visualised for which = 1: this is based on an auxiliary run of Gibbs sampling for a fixed model")
#     if(show[2] & (!x$G.sel)) warning("optimal model visualised for which = 2 ")
#     if(show[2]){ 
#       x$Z <- x$Z[[which.max(table(x$samples$G))]]
#       x$count <- rep(1, nrow(x$Z))
#       if(is.null(x$classprob)) x$classprob <- colMeans(x$Z)
#       }
#     if(any(show[3:4])){ 
#       warning("density estimates not supported for collapsed Gibbs sampling algorithm")
#       show[3:4]<- FALSE
#     }
#     which<- c(1:4,10)[show]
#     NextMethod("plot")
#   }
