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
			
	Last modification of this code: Mon 09 May 2016 15:38:03 IST   */

#ifndef __BLCA_INITIALIZE_H__
#define __BLCA_INITIALIZE_H__

#include "BLCA_required_libs.h"
#include "BLCA_mixmod.h"

void set_prior_on_number_of_components(struct mix_mod *mixmod,int type) ;

int initialize_simple(struct mix_mod *mixmod,int numgroups) ;

#endif
