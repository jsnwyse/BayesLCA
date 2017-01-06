plot.blca.collapsed <-
  function(x, which=1L, main="", ...){
    #class(x)<- "blca"
    #print("NextMethodUsed")
    show<- rep(FALSE, 5)
    show[which]<- TRUE
    if(show[1] & is.null(x$itemprob)){
      warning("Item probabilities not calculated. Use ?blca.collapsed.post.hoc.estimates for more information.")
      show[1]<- FALSE 
    }
    if(show[1] & (!x$fixed.G | x$variable.selection)) warning("Optimal model only is visualised for which = 1.")
    if(show[2] & (!x$fixed.G)) warning("Optimal model only is visualised for which = 2.")
    if(show[2]){ 
      x$Z <- x$Z[[which.max(table(x$samples$Giter))]]
      x$count <- rep(1, nrow(x$Z))
      if(is.null(x$classprob)) x$classprob <- colMeans(x$Z)
      }
    if(any(show[3:4])){ 
      warning("Density estimates not supported for collapsed Gibbs sampling algorithm.")
      show[3:4]<- FALSE
    }
    which<- c(1:4,10)[show]
    NextMethod("plot")
  }
