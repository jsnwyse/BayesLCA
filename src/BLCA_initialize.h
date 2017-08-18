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

void BLCA_set_prior_on_number_of_components(struct mix_mod *mixmod,int type) ;

void BLCA_initialize_data( struct mix_mod *mixmod, int *Y );

void BLCA_initialize_priors( struct mix_mod *mixmod, double *alpha_prior, double *beta_prior, int type );

void BLCA_initialize_Gibbs_sampler( int init_type, struct mix_mod *mixmod ) ;

void BLCA_initialize_EM( int init_type, double *init_vals, struct mix_mod *mixmod , double *group_weights, double *prob_variables) ;

void BLCA_initialize_VB( int init_type, double *init_vals, struct mix_mod *mixmod , double *alpha_ud, double *beta_ud );

void BLCA_initialize_Boot( struct mix_mod *mixmod, double *group_weights, double *prob_variables  );

#endif
