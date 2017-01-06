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
	
#ifndef __BLCA_VARIABLE_SELECTION_H__
#define __BLCA_VARIABLE_SELECTION_H__

#include "BLCA_mixmod.h"

void BLCA_update_model_by_variable_include_exclude(struct mix_mod *mixmod,int *accepted,int *proposed,int var) ;

#endif
