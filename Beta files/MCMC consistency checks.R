library(BayesLCA)

## check functionality MCMC functions

data(Alzheimer)

alz <- data.blca(Alzheimer)

##-----------------
## blca.gibbs - old
##-----------------

fit <- blca.gibbs(Alzheimer, G = 2)
fit <- blca.gibbs(alz, G = 2)

G = 2

fit <- blca.gibbs(alz, G = 2, alpha = rep(1, G), beta = rep(1, G), delta = 0.5, start.vals = "prior", thin = 0.3, burn.in = 500, iter = 15500, verbose.update = 1000)

##-----------------
## collapsed Gibbs
##-----------------


fit <- blca.collapsed(Alzheimer, G = 2)
fit <- blca.collapsed(alz, G = 2) ## ERROR due to data.blca format

G = 2

## no start vals, or verbose arguments  
## thin is other way around i.e, thin = 3, not thin = 1/3

fit <- blca.collapsed(Alzheimer, G = 2, alpha = rep(1, G), beta = rep(1, G), delta = 0.5, thin = 3, burn.in = 500, iter = 15500)

fit <- blca.collapsed(Alzheimer, G = 2, alpha = rep(1, G), beta = rep(1, G), delta = 0.5, thin = 3, burn.in = 500, iter = 15500, post.hoc.est = FALSE, variable.selection = TRUE)

