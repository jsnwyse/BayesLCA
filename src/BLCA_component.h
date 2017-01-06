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
	
#ifndef __BLCA_COMPONENT_H__
#define __BLCA_COMPONENT_H__

#include "BLCA_required_libs.h"
#include "BLCA_mixmod.h"

void BLCA_allocate_component(struct component *component,struct mix_mod *mixmod) ;

void BLCA_free_component(struct component *component,struct mix_mod *mixmod) ;

void BLCA_copy_component(struct component *component_original,struct component *component_target,struct mix_mod *mixmod) ;

void BLCA_recompute_sufficient_statistics_for_components(struct mix_mod *mixmod) ;

double BLCA_compute_log_data_probability_with_inclusion_in_component(int *x,struct component *component,struct mix_mod *mixmod) ;

double BLCA_compute_log_marginal_likelihood_with_inclusion_in_component(int *x,struct component *component,struct mix_mod *mixmod) ;

double BLCA_compute_log_data_probability_component(struct component *component,struct mix_mod *mixmod) ;

void BLCA_recompute_marginal_likelihood_component(struct component *component,struct mix_mod *mixmod) ;

void BLCA_recompute_marginal_likelihood_undiscriminating_variables(struct component *undiscriminating,struct mix_mod *mixmod) ;

#endif
