plot.blca.multicat <-
  function(x, which=1L, main="", col1=heat.colors(12), by1 = c("variable", "group"), auto.layout = TRUE, ...){
    logistore<- devAskNewPage()
    devAskNewPage(TRUE)
    
    if( any(which>2) ){
    	warning("\n Plotting for more than two categories: item probability and classification uncertainty plot only.")
    } 
    
    show<- rep(FALSE, 11)
    show[which]<-TRUE
    #print(show)
if(show[1]){
  
  #by1 <- "variable"
  #by1 <- "group"
  by1 <- match.arg(by1)
  old_mfrow <- par()$mfrow
  
  z<-NULL
  
  btau<- x$classprob
  btheta<- x$itemprob.tidy
  
  G <- length(btau)
  names(btau)<- paste("Group", 1:G)
  
  z$ItemProb<- btheta
  z$MembershipProb<- btau
  
  M <- nlevels(btheta$variable)
  lev.var <- levels(btheta$variable)
  
  olevel <- apply(t(1:M), 2, function(x) min(which(btheta$variable == lev.var[x])))
  btheta$variable <- factor(btheta$variable, levels = lev.var[order(olevel)])
  lev.var <- levels(btheta$variable)
  
  if(by1 == "variable"){
    if (auto.layout) {
      mfrow <- coda:::set.mfrow(Nchains = 1, Nparms = M, nplots = 1)
      oldpar <- par(mfrow = mfrow)
    }
    for(m1 in 1:M){
      #cat("\n", paste(lev.var[m1]), ":\n")
      temp1 <- btheta[btheta$variable == lev.var[m1], ]
      temp1$category <- factor(temp1$category)
      nlev.temp <- nlevels(temp1$category)
      temp1 <- matrix(temp1[, 1], nrow = G, ncol = nlev.temp,  dimnames = list(names(btau), levels(temp1$category)))
      
      #C1 <- ncol(temp1)#ncol(Theta)
      #G<- nrow(Theta)
      #xcut<-c(0,cumsum(btau))
      #ycut<-0:C1
      
      #image.plot(ycut, xcut, t(temp1), axes=FALSE, zlim=c(0,1), ylab="Groups", xlab="Categories", main=main, col=col1)
      mosaicplot(temp1, color = 1:nlev.temp + 1, ylab="Categories", xlab="Groups", main = main)
      
      mtext(lev.var [m1], 3, 0.25)
    }
  } else{
    if(by1 == "group"){
      
      dum.mat <- matrix(0, nrow = M, ncol = nlevels(btheta$category))
      dimnames(dum.mat) <- list(lev.var, levels(btheta$category))
      if (auto.layout) {
        mfrow <- coda:::set.mfrow(Nchains = 1, Nparms = G, nplots = 1)
        oldpar <- par(mfrow = mfrow)
      }
      for(g1 in 1:G){
        #cat("\n", paste(names(btau)[g1]), ":\n")
        temp1 <- btheta[btheta$group == names(btau)[g1], ]
        for(m1 in 1:M){
          sel1 <- temp1$variable == lev.var[m1]
          dum.mat[m1, 1:length(temp1$category[sel1])] <- temp1[sel1, 1]
        }
        
        mosaicplot(dum.mat, color = 1:nlevels(btheta$category)+1, xlab = "Variables", ylab = "Item Probabilities", main = main)
        mtext(paste("Group", g1), 3, 0.25)
        
      } 
    }
  }
}#Show2 Multicat Parameter Mosaicplot
    old_mfrow <- par()$mfrow
    show[1]<- FALSE
    which <- which(show)
    plot.blca(x, which, main = main, col1 = col1, ...)
  }
