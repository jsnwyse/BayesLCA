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


#include "BLCA_component.h"

void BLCA_allocate_component(struct component *component,struct mix_mod *mixmod)
{
	/*allocate memory for a component*/
	int j;
	
	component->N = calloc(mixmod->d,sizeof(int *));
	
	for(j=0;j<mixmod->d;j++){
		component->N[j] = calloc(mixmod->ncat[j],sizeof(int));
	}
	
	if( !mixmod->collapsed || mixmod->EM_fit )
	{
		component->prob_variables = calloc(mixmod->d,sizeof(double *));
		for(j=0;j<mixmod->d;j++){
			component->prob_variables[j] = calloc(mixmod->ncat[j],sizeof(double));
		}
	}
	
	return;
}

void BLCA_free_component(struct component *component,struct mix_mod *mixmod)
{
	
	int j;
	
	for(j=0;j<mixmod->d;j++){
		free(component->N[j]);
	}
	free(component->N);
	
	if( !mixmod->collapsed || mixmod->EM_fit )
	{
		for(j=0;j<mixmod->d;j++){
			free(component->prob_variables[j]);
		}
		free(component->prob_variables);
	}
	
	return;
	
}

void BLCA_copy_component(struct component *component_original,struct component *component_target,struct mix_mod *mixmod)
/*copy the contents of the first argument into the second component argument*/
{

	int i,j;
	
	component_target->n_g = component_original->n_g;
	for(i=0;i<mixmod->d;i++){
		for(j=0;j<mixmod->ncat[i];j++){
			component_target->N[i][j] = component_original->N[i][j];
		}
	}

	component_target->log_prob = component_original->log_prob;
	
	return;

}

void BLCA_recompute_sufficient_statistics_for_components(struct mix_mod *mixmod)
{

	int i,j,k;
	
	/*undiscriminating variables*/
	
	
	
	/*discriminating variables*/
	
	for(k=0;k<mixmod->G;k++){	
		/*reset*/
		for(j=0;j<mixmod->d;j++){
			for(i=0;i<mixmod->ncat[j];i++){
				mixmod->components[k]->N[j][i] = 0;
				mixmod->undiscriminating->N[j][i] = 0;
			}
		}
		mixmod->components[k]->n_g = 0;
	}
	
	int z;
	
	for(i=0;i<mixmod->n;i++){
	
		z = mixmod->z[i];	
		
		mixmod->components[ z ]->n_g += 1;
		
		for(j=0;j<mixmod->d;j++){
			mixmod->components[ z ]->N[j][ mixmod->Y[j][i] ] += 1;
			mixmod->undiscriminating->N[j][ mixmod->Y[j][i] ] += 1;
		}
	
	
	}
	
	return;
}


double BLCA_compute_log_data_probability_with_inclusion_in_component(int *x,struct component *component,struct mix_mod *mixmod)
{

	int i,j,I;
	double log_prob = 0.;
	
	for(j=0;j<mixmod->d;j++){
	
		if(mixmod->varindicator[j]){
	
			log_prob += lgamma(mixmod->ncat[j] * mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) - lgamma( component->n_g+1. + mixmod->ncat[j]*mixmod->beta );
		
			for(i=0;i<mixmod->ncat[j];i++){
		
				I = (x[j] == i) ?  1 : 0 ;
		
				log_prob += lgamma(component->N[j][i] + I + mixmod->beta);
		
			}
		
		}
		
	}

	/*if(isnan(log_prob)){
		mixmod_warning(1);
		printf("\nFrom compute_log_data_probability_with_inclusion_in_component. Result = %.10f,n_g = %d\n,sum_sq_norm = %.10f,sq_norm = %.10f",sum_sq_norm - sq_norm/(component->n_g+kappa) + kappa*xi2 + gamma,component->n_g,sum_sq_norm,sq_norm);
	}*/

	return(log_prob);
}

double BLCA_compute_log_marginal_likelihood_with_inclusion_in_component(int *x,struct component *component,struct mix_mod *mixmod)
{
	/*function used for class prediction*/

	int i,j,I;
	double log_prob = 0.;
	
	log_prob = lgamma(component->n_g+1.+mixmod->alpha);
	
	for(j=0;j<mixmod->d;j++){
	
		if(mixmod->varindicator[j]){
	
			log_prob += lgamma(mixmod->ncat[j] * mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) - lgamma( component->n_g + 1. + mixmod->ncat[j]*mixmod->beta );
		
			for(i=0;i<mixmod->ncat[j];i++){
		
				I = (x[j] == i) ?  1 : 0 ;
		
				log_prob += lgamma(component->N[j][i] + I + mixmod->beta);
		
			}
		
		}
		
	}

	return(log_prob);
}





double BLCA_compute_log_data_probability_component(struct component *component,struct mix_mod *mixmod)
{

	int i,j;
	double log_prob = 0.;

	for(j=0;j<mixmod->d;j++){
	
		if(mixmod->varindicator[j]){
		
			log_prob += lgamma(mixmod->ncat[j] * mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) - lgamma(component->n_g + mixmod->ncat[j]*mixmod->beta );
		
			for(i=0;i<mixmod->ncat[j];i++){
		
				log_prob += lgamma(component->N[j][i] + mixmod->beta);
		
			}
		
		}
		
	}
	
	/*if(isnan(log_prob)){
		mixmod_warning(1);
		printf("\nFrom compute_log_data_probability_component\n");
	}*/

	return(log_prob);	


}


void BLCA_recompute_marginal_likelihood_component(struct component *component,struct mix_mod *mixmod)
{
	/*computing these quantities will save some time for the Gibbs update of the labels*/
	int i,j;
	double log_prob;
	
	log_prob = lgamma(component->n_g+mixmod->alpha);
	
	if(component->n_g > 0){
	
		for(j=0;j<mixmod->d;j++){
		
			if(mixmod->varindicator[j]){
			
				log_prob += lgamma(mixmod->ncat[j] *mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) - lgamma(component->n_g + mixmod->ncat[j]*mixmod->beta);
			
				for(i=0;i<mixmod->ncat[j];i++){
			
					log_prob += lgamma(component->N[j][i] + mixmod->beta);
			
				}
			
			}
			
		}
		
	}
	
	/*if(isnan(log_prob)){
		mixmod_warning(1);
		printf("\nFrom recompute_marginal_likelihood_component\n");
	}*/
	
	component->log_prob = log_prob;

	return;
}


void BLCA_recompute_marginal_likelihood_undiscriminating_variables(struct component *undiscriminating,struct mix_mod *mixmod)
{

	int i,j;
	double log_prob;
	
	log_prob = 0.;
	
	for(j=0;j<mixmod->d;j++){
	
		if(!mixmod->varindicator[j]){ /*variables not in mixture model*/
		
			log_prob += lgamma(mixmod->ncat[j]*mixmod->beta) - mixmod->ncat[j]*lgamma(mixmod->beta) - lgamma(undiscriminating->n_g + mixmod->ncat[j]*mixmod->beta);
			
			for(i=0;i<mixmod->ncat[j];i++){
			
				log_prob += lgamma(undiscriminating->N[j][i] + mixmod->beta);
			
			}
		
		}
		
	}
	
	undiscriminating->log_prob = log_prob;
	
	return;
	
}




