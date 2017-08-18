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

#ifndef __BLCA_RESULTS_H__
#define __BLCA_RESULTS_H__

#include "BLCA_utils.h"

struct results
/*a structure to store all of the results of the analysis*/
{
	
	int *ngroups; /*gives number of groups at each iteration*/
	
	int **memberships; /*gives membership at each iteration*/
	
	int *MAP_memberships; /*gives membership with maximum prob for each i=1,\dots,n*/
	
	int **variable_indicator; /*gives binary vector of variable inclusion for each iteration*/
	
	double *variable_prob_inclusion;
	
	double *log_posterior; /*stores the log posterior for all iterations...*/
	
	int nsteps;
	
	double tol_obtained;
	
	int proposed_m1;
	
	int accepted_m1;
	
	int proposed_m2;
	
	int accepted_m2;
	
	int proposed_m3;
	
	int accepted_m3;
	
	int proposed_eject;
	
	int accepted_eject;
	
	int proposed_absorb;
	
	int accepted_absorb;
	
	int proposed_remove_variable;
	
	int accepted_remove_variable;
	
	int proposed_add_variable;
	
	int accepted_add_variable;
	
	int proposed_include_exclude_variable;
	
	int accepted_include_exclude_variable;
	
};

void BLCA_allocate_results(struct results *results,int iterations,int burn_in,int thin_by,int len,int d) ;

void BLCA_reset_results(struct results *results);

void BLCA_free_results(struct results *results,int iterations,int burn_in,int thin_by) ;

void BLCA_allocate_results_x2(struct results *results,int iterations,int burn_in,int thin_by,int len,int d,int writeToFile) ;

void BLCA_free_results_x2(struct results *results,int iterations,int burn_in,int thin_by,int writeToFile) ;

int BLCA_write_out_results(struct results *results,int N,int datasize,int datadimension,int onlyGibbs,int fixedG,int selectVariables) ;

#endif
