/*Functions for the fitting of Latent Class Analysis models using MCMC
	methods. Two implementations of the model are included in the Bayesian
	formulation: collapsed and not collapsed.
	
	Author:	Jason Wyse,
			School of Computer Science and Statistics,
			Lloyd Institute,
			Trinity College,
			Dublin 2,
			Ireland.
			mailto: wyseja@tcd.ie
			
	Last modification of this code: Thu 04 Aug 2016 21:43:55 IST    */
	
#ifndef __BLCA_DENSITY_H__
#define __BLCA_DENSITY_H__

#include "BLCA_mixmod.h"

double BLCA_log_normalizing_constant_model(int G,struct mix_mod *mixmod) ;

double BLCA_l_prior_variable_include(int D,struct mix_mod *mixmod) ;

double BLCA_get_full_log_posterior(struct mix_mod *mixmod) ;

double BLCA_get_full_log_posterior_x2(struct mix_mod *mixmod) ;

double BLCA_get_log_likelihood(struct mix_mod *mixmod) ;

double BLCA_get_VB_bound( struct mix_mod *mixmod ) ;

#endif
