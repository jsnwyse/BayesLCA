library(BayesLCA)

## check functionality MCMC functions

set.seed( 1984 )

##---------------------
##--- Alzheimer data---
##---------------------

data(Alzheimer)

alz <- data.blca(Alzheimer)

##-----------------
## blca.gibbs - old
##-----------------

fit <- blca.gibbs(Alzheimer, G = 2)
fit <- blca.gibbs(alz, G = 2)

G <- 2

fit <- blca.gibbs(alz, G = 2, alpha = rep(1, G), beta = rep(1, G), delta = 0.5, start.vals = "prior", thin = 0.3, burn.in = 500, iter = 15500, verbose.update = 1000)

##-----------------
## collapsed Gibbs
##-----------------


fit <- blca.collapsed(Alzheimer, G = 2, G.sel=T)
fit <- blca.collapsed(alz, G = 2) 

G = 2


fit <- blca.collapsed(Alzheimer, G = 2, alpha = rep(1, G), beta = rep(1, G), delta = 0.5, thin = .3, burn.in = 500, iter = 15500, G.sel=T )


fit <- blca.collapsed(Alzheimer, G = 2, G.max=10, alpha = rep(1, G.max), beta = rep(1, G.max), delta = 0.5, thin = .3, burn.in = 500, iter = 15500, G.sel=T, var.sel=T )



## ----------------------
## --- simulated data ---
## ----------------------

## Example 2  Dean and Raftery (2010)


#this is a three class model
ncat = c(3,2,4,3,3,4,5,2,3,4)

#matrix of probabilities

Prob = c(.1,.3,.6,
         .1,.5,.2,
         .8,.2,.2,
         .5,.1,.7,
         .5,.9,.3,
         .2,.7,.2,
         .2,.1,.6,
         .3,.1,.1,
         .3,.1,.1,
         .1,.6,.4,
         .5,.1,.4,
         .4,.3,.2,
         rep(.4,3),
         rep(.5,3),
         rep(.1,3),
         rep(.2,3),
         rep(.4,3),
         rep(.1,3),
         rep(.3,3),
         rep(.2,3),
         rep(.3,3),
         rep(.3,3),
         rep(.1,3),
         rep(.1,3),
         rep(.2,3),
         rep(.8,3),
         rep(.7,3),
         rep(.1,3),
         rep(.2,3),
         rep(.1,3),
         rep(.2,3),
         rep(.1,3),
         rep(.6,3))

P = matrix(Prob,nrow=sum(ncat),ncol=3,byrow=TRUE)

itpr <- list()
for( j in 1:length(ncat) )
{
	if( j==1 ) idx <- 0 else idx <- sum(ncat[1:(j-1)])
	itpr[[j]] <- t( P[ (idx+1):(idx+ncat[j]), ] )
}

w = c(.3,.4,.3)

Y <- rlca( n=1000, classprob=w, itemprob=itpr, ncat=ncat )


## Gibbs sampler ###
fit <- blca.gibbs( Y, ncat=ncat, G = 3)

## EM algorithm for fitting ###
fit <- blca.em( Y, ncat=ncat, G = 3)

## VB algorithm ###
fit <- blca.vb( Y, ncat=ncat, G = 3)








