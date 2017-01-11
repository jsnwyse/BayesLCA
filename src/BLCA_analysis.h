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
			
	Last modification of this code: Thu 04 Aug 2016 21:02:24 IST    */

#ifndef __BLCA_ANALYSIS_H__
#define __BLCA_ANALYSIS_H__

#include "BLCA_mixmod.h"

struct results *BLCA_analysis_MCMC_collapsed( struct mix_mod *mixmod, int num_iteration, int num_burnin, int thin_by, int fixedG, int onlyGibbs, 
				int selectVariables, int *group_memberships, int *variable_inclusion_indicator, int *n_groups, double *log_joint_posterior,
				 double *prior_include, int *var_pattern );

void BLCA_analysis_MCMC_Gibbs_sampler( struct mix_mod *mixmod, int num_iteration, int num_burnin, int thin_by, int *group_memberships, 
				double *group_weights, double *variable_probabilities, double *log_joint_posterior, int verbose, int verbose_update );

void BLCA_analysis_EM( struct mix_mod *mixmod, int max_iterations, int *iterations, double *membership_probabilities, double *group_weights, 
			double *variable_probabilities, double *log_likelihood, int MAP, double tol, int *converged ) ;

void BLCA_E_step( struct mix_mod *mixmod );

void BLCA_M_step( struct mix_mod *mixmod );

void BLCA_analysis_VB( struct mix_mod *mixmod, int max_iterations, int *iterations, double *phi, double *alpha_ud, 
				double *beta_ud, double *lower_bound, double tol, int *converged ) ;

void BLCA_VB_phi_step( struct mix_mod *mixmod );

void BLCA_VB_alpha_beta_step( struct mix_mod *mixmod );

#endif
