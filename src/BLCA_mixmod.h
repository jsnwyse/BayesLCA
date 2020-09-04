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

#ifndef __BLCA_MIXMOD_H__
#define __BLCA_MIXMOD_H__

#include "BLCA_required_libs.h"
#include "BLCA_utils.h"
#include "BLCA_digamma_func.h"

struct component
{
	int in_use; /*indicator saying whether the component is in use or not*/
	
	int n_g; /*number of members of component*/
	
	//int *indexes_g; /*indexes of members of component in Y (row-wise)*/
	
	int **N; /*gives the number in category k of variable j*/
	
	double **prob_variables; /*matrix of variable probabilities*/
	
	double **beta_ud; /*dirichlet parameters for VB*/
	double **di_beta_ud;
	double *di_sum_beta_ud;
	
	double log_prob; /*quantity that can be updated to save computations*/

};

struct mix_mod
/*structure to hold mixture essentials*/
{
	int collapsed; /*logical indicating collapsed model or not*/

	int G; /*total number of groups*/
	
	int n; /*total number of data*/
	
	int n_gibbs; /*number of data pts to update in Gibbs*/
	
	int d; /*total number of variables recorded*/
	
	int maxgroups;
	
	int *ncat; /*number of categories for each of the variables*/
	
	int *varindicator; /*indicates whether the variable is in or not*/

  int *y; /*pointer to data passed from R*/
	
	int **Y; /*the raw data stacked as variable \times obs (J \times N)*/
	
	int **Yobs; /*a second copy of the raw data stacked as obs \times variable */
	
	/*two copies of the data are stored as I think its better for memory access
		- we can use one or the other depending on the calculation to be carried out*/
	
	int *z; /*group memberhsip indicator*/
	
	double **s; /*group probability membership*/
	
	int *whereis; /*stores the index in components of a component*/

	struct component **components; /*pointer to array of pointers to components*/
	
	struct component *undiscriminating; /*for variable selection: pointer to component for 
														undiscriminating variables*/
	
	/*other hyperparameters*/
	
	int prior_type; //this will tell us whether prior is identifiable or not- which computation to do
	
	double *alpha_prior; // vector of prior on weights
	
	double lg_alpha_sum; // sum of gamma
	
	double lg_sum_alpha;
	
	double ***beta_prior; // array (possibly varying in size) of Dirichlet prior on weights 
	
	double **lg_beta_sum; // array with sum of log gamma of betas
	
	double **lg_sum_beta; // array with log gamma of sum of betas
	
	double alpha; /*alpha: dirichlet prior on weights symmetric*/
	
	double beta; /*beta: dirichlet prior on the category probabilities symmetric*/
	
	double prior_prob_variable_include; /*prior probability of inclusion of any of the variables*/ 
	
	double *log_prior_G; /*this is the prior for the number of groups... poisson(1)*/
							
	double *table_a; /*this is a lookup table for the values of a when ejecting or
								absorbing components*/
	
	double *weights; /*component weights if collapsed = FALSE*/
	
	double *alpha_ud; /*dirichlet parameters for the variational approximation*/
	double *di_alpha_ud;
	double di_sum_alpha_ud;
	
	int hprior_model; /*logical: sample a hyperprior on the probability of variable inclusion*/
	
	double hprior_model_a0; /*hyperparameter values for the Beta hyperprior if hprior_model == TRUE*/
	
	double hprior_model_b0;
	
	double log_like;
	
	double log_prior;
	
	FILE *fp_log; /*file pointer for the debugger log*/
	
	int EM_fit ;
	
	int EM_MAP; /*find the map in an EM fit?*/ 
	
	//int EM_SMALL_PROB_WARNING;
	
	int VB;
	
	int BOOT;
	
	int *boot_idx;
	
	int component_compute; // the component currently being computed

};

#include "BLCA_component.h"
#include "BLCA_initialize.h"
#include "BLCA_results.h"
#include "BLCA_post_hoc.h"
#include "BLCA_label_updates.h"
#include "BLCA_eject_absorb.h"
#include "BLCA_variable_selection.h"
#include "BLCA_density.h"
#include "BLCA_analysis.h"


struct mix_mod *BLCA_allocate_mixmod(int n, int d, int G_max, int G, double *prior_hparams,int *ncat, 
                                     int COLLAPSED, int EM_fit, int EM_MAP, int VB, int BOOT );

struct mix_mod *BLCA_clone_mixmod( struct mix_mod *mm );

void BLCA_free_mixmod(struct mix_mod *mm);

#endif
